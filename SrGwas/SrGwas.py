from miscSupports import validate_path, directory_iterator, load_yaml, FileOut, terminal_time, chunk_list, flatten
from csvObject import CsvObject, write_csv
from bgen_reader import custom_meta_path
from pysnptools.distreader import Bgen
import statsmodels.api as sm
from random import sample
from pathlib import Path
from scipy import stats
import pandas as pd
import numpy as np
import re


class SrGwas:
    def __init__(self, args):

        # Load the args from the yaml file
        self.args = load_yaml(args)
        self.write_dir = self.args["output_directory"]

        # Set the gen file info, set the output path for the memory files, and load the file reference
        self.gen_type = self.args["gen_type"]
        self.gen_directory = self.args["path_to_gen_files"]
        self.target_chromosome = self.args["target_chromosome"]
        self.file_name = f"{self.args['output_name']}_Chr{self.target_chromosome}"
        custom_meta_path(validate_path(self.args["memory_file_location"]))

        # Setup logger
        self.logger = FileOut(self.write_dir, self.file_name, "log", True)
        self.logger.write(f"Setup {terminal_time()}")

        # Variable info, load the genetic reference, and sort both it and the external variables so they match on iid
        self.phenotype = self.args["phenotype"]
        self.covariant = self.args["covariant"]
        self.gen, self.df, self.genetic_iid = self._setup_variables()
        self.total_obs = len(self.df)
        self.logger.write(f"Set {self.gen.iid_count} in Genetic file and {len(self.df)} in variable file for "
                          f"{self.phenotype}~{self.covariant}")

        # Check that we only have a single version of phenotypic columns, if the file contained one of these names this
        # could by why we now have duplicates
        if len(self.df[f"{self.phenotype}RES"].shape) > 1:
            self.logger.write(f"Found a duplicated column for phenotypic residuals, removing")
            self.df = self.df.loc[:, ~self.df.columns.duplicated()]

        # Set output file
        if self.args["method"] != "set_snp_ids":
            self.output = FileOut(validate_path(self.write_dir), self.file_name, "csv")
            headers = [[f"M{i}_{h}" for h in ["coef", "std_err", "pvalue", "obs", "r2", "chi2tail", "95%lower",
                                              "95%upper"]]
                       for i in range(1, 5)]
            self.output.write_from_list(["Snp"] + flatten(headers))

        # System args
        self.iter_size = self.args["array_size"]
        self.start_index = 0

        # Start the method that has been assigned if method has been set
        if self.args["method"]:
            getattr(self, self.args["method"])()

        self.logger.write(f"Finished predefined {terminal_time()}")

    def __repr__(self):
        return f"SrGwas object Controller"

    def _setup_variables(self):
        """
        The order of IID in genetic file may not equal to submission, this sorts the arrays to be equivalent.

        :return: Bgenfile for this chromosome as well as a pandas dataframe of the external variables
        """

        # Load the variables as pandas dataframe and setup the reference genetic file for this chromosome
        df = pd.read_csv(validate_path(self.args["variables"]))
        gen = Bgen(self._select_file_on_chromosome())
        self.logger.write(f"...Loaded external variables {terminal_time()}")

        # Validate that the variables we have set in the formula exist in the DataFrame
        [self._validate_variable(df, cont, "Continuous") for cont in self.covariant]
        assert self.args["phenotype"], "GWAS requires a phenotype"

        # Recast IID as an int
        df["IID"] = [self._strip_iid(iid) for iid in df["IID"].tolist()]

        # Isolate the IID to match against the variables IID and create the reference
        genetic_iid = np.array([self._strip_iid(iid) for _, iid in gen.iid])
        genetic_position = gen.iid

        # Remove any IID that is in the external data array but not in the genetic array
        out = np.in1d(df["IID"].to_numpy(), genetic_iid)
        df = df[out]

        # Remove any IID that is in the genetic array but not in the external data
        out = np.in1d(genetic_iid, df["IID"].to_numpy())
        genetic_iid = genetic_iid[out]
        genetic_position = genetic_position[out]

        # Sort both arrays to be in the same order
        df = df.sort_values(by=['IID'], ascending=True)
        gen = gen[gen.iid_to_index(genetic_position[np.argsort(genetic_iid)]), :]

        for index, v in enumerate(df.columns):
            if v in [self.phenotype] + [self.covariant]:
                df[v] = df[v].apply(pd.to_numeric)

        # Create an IID array of the genetic iid
        genetic_iid = pd.DataFrame(genetic_iid)
        genetic_iid.columns = ["IID"]

        # Add a constant and the residualised phenotype to the databases
        df["Constant"] = [1 for _ in range(len(df))]
        self.covariant = self.covariant + ["Constant"]
        result = sm.OLS(df[self.phenotype], df[self.covariant], missing='drop').fit()
        df = pd.concat([df, pd.DataFrame(result.resid, columns=[f"{self.phenotype}RES"])], axis=1)

        # Remove non used data to save memory
        return gen, df[["IID", self.phenotype, f"{self.phenotype}RES"] + self.covariant + ["Constant"]], genetic_iid

    def _select_file_on_chromosome(self):
        """
        For a given chromosome, get the respective file from the genetic directory

        :return: Path to the current file as a string representation of a Path from pathlib
        :rtype: str

        :raises IndexError: If not file is found
        """
        for file in directory_iterator(self.gen_directory):
            if Path(self.gen_directory, file).suffix == self.gen_type:
                try:
                    if int(re.sub(r'[\D]', "", Path(self.gen_directory, file).stem)) == self.target_chromosome:
                        return str(Path(self.gen_directory, file).absolute())
                except (ValueError, TypeError):
                    continue

        raise IndexError(f"Failed to find any relevant file for {self.target_chromosome} in {self.gen_directory}")

    @staticmethod
    def _validate_variable(variables, v, var_type):
        """Check the variable exists within the columns"""
        if v and v != "null":
            assert v in variables.columns, f"{var_type} variable {v} not in variables {variables.columns}"

    @staticmethod
    def _strip_iid(iid):
        """Strip IID of any non numeric characters"""
        return int(re.sub(r'[\D]', "", str(iid)))

    def _select_snps(self):
        """
        We may only want to run a subset of snps. If so, then this loads the snp indexes from a csv. Else, just return
        all the snp ids

        :return: A list of snp ids
        :rtype: list[snp]
        """
        if self.args["snps"]:
            return CsvObject(validate_path(self.args["snps"]), set_columns=True, column_types=int)[0]
        else:
            return [i for i in range(self.gen.sid_count)]

    def set_snp_ids(self):
        """
        Isolate a subset of snps based on pre-defined named snps in a csv, passed as a str to snps_to_id, or a random
        set of snps of total == pre-defined int, where the int is set to snps_to_id.

        :return: Nothing, write the id's to a csv then stop
        :rtype: None
        :raise TypeError: If a str/ int is not passed
        """

        if isinstance(self.args["snps_to_id"], str):
            snps = CsvObject(validate_path(self.args["snps_to_id"]), set_columns=True)[0]
            snp_list = [f"{snp},{snp}" for snp in snps]
            snp_list = self.gen.sid_to_index(snp_list).tolist()

        elif isinstance(self.args["snps_to_id"], int):
            snp_list = sample(range(self.gen.sid_count), self.args["snps_to_id"])

        else:
            raise TypeError(f"set_snp_ids expects an int or a str yet was passed {type(self.args['snps_to_id'])}\n"
                            f"If you want to run a GWAS with all snps, simply set snps: null and ignore this method")

        write_csv(self.write_dir, f"{self.file_name}_Snps", ["Snp"], snp_list)
        self.logger.write(f"Constructed snp id list {terminal_time()}")

    def gwas(self):
        """
        Create genetic residuals by regressing your covariant on the snp or run a more traditional gwas of

        phenotype ~ dosage + covariant_1 + ... covariant_N

        :return: Nothing, write line to fine when residuals have been estimated
        :rtype: None
        """
        # Isolate which snps are to be used
        snp_ids = self._select_snps()
        snp_chunk_list = chunk_list(snp_ids[self.start_index:], self.iter_size)

        for chunk_id, snp_chunk in enumerate(snp_chunk_list, 1):
            self.logger.write(f"Chunk {chunk_id} of {len(snp_chunk_list)}")

            # Instance the memory for all individuals (:) for snp i
            current_snps = self.gen[:, snp_chunk]

            # Transform bgen dosage of [0, 1, 0] -> 0, 1, or 2 respectively.
            dosage = sum(np.array([snp * i for i, snp in enumerate(current_snps.read(dtype=np.int8).val.T)],
                                  dtype=np.int8))
            self.logger.write(f"Loaded Chunk {chunk_id}: {terminal_time()}")

            # Isolate the snp names
            snp_names = [snp.split(",")[1] for snp in current_snps.sid]

            # Construct a dataframe from the demeaned covariant dataframe and these snps
            snp_df = pd.DataFrame(dosage).T
            snp_df.columns = snp_names

            # Merge this snp data on IID
            snp_df = pd.concat([self.genetic_iid, snp_df], axis=1)
            df = self.df.merge(snp_df, left_on="IID", right_on="IID")

            for i, snp in enumerate(snp_names):
                if i % (self.iter_size / 10) == 0:
                    self.logger.write(f"snp {i}/{len(snp_chunk)}: {terminal_time()}")

                # Define the output list
                out_list = [snp]

                # Model 1: Traditional OLS
                result = sm.OLS(df[self.phenotype], df[[snp] + self.covariant], missing='drop').fit()
                out_list = out_list + self.results_out(result, snp, len(self.covariant) + 1)

                # Model 2: Phenotypic Residual
                result = sm.OLS(df[f"{self.phenotype}RES"], df[[snp, "Constant"]], missing='drop').fit()
                out_list = out_list + self.results_out(result, snp, 2)

                # Model 3: Genetic residual
                g_res = sm.OLS(df[snp], df[self.covariant], missing='drop').fit()
                g_res = pd.concat([pd.DataFrame(g_res.resid, columns=[snp]), df["Constant"]], axis=1)

                result = sm.OLS(df[self.phenotype], g_res, missing='drop').fit()
                out_list = out_list + self.results_out(result, snp, 2)

                # Model 4: Genetic residual on phenotypic residuals
                result = sm.OLS(df[f"{self.phenotype}RES"], g_res, missing='drop').fit()
                out_list = out_list + self.results_out(result, snp, 2)
                self.output.write_from_list(out_list, True)

    def results_out(self, results, v_name, model_k, alpha=0.05):
        """
        Returns for each variable in the list of variables

        [Parameters, standard error, p values, obs, 95%min CI, 95%max CI]

        Notes
        -----
        To make models cross comparable results are adjusted

        For the coefficients, we use the squared T statistics, using the sample size N as the denominator in the
        variance estimator (RSS/N) instead of RSS/(N-k_2). Standard errors are adjusted based the number of covariants
        selected, rather than the n - k -1 of the model. Confidence intervals are adjusted on the standard normal
        distribution.

        :param results: The mostly unadjusted results from OLS bar the degrees of freedom that was adjusted for clusters
        :type results: statsmodels.regression.linear_model.RegressionResults

        :param v_name: A string of the variable to extract from
        :type v_name: str

        :param model_k: The n-k of this model to adjust the standard errors
        :type model_k: int

        :param alpha: The Significance level for the confidence interval, defaults to 0.05 for a 95% confidence interval
        :type alpha: float

        :return: A list of lists, where each list are the results in float
        :rtype: list[list[float, float, float, float, float, float]]
        """

        # Adjust the coefficient
        snp_estimate = results.params[v_name]
        snp_variance = results.cov_params()[v_name][v_name]
        estimate_adj = (snp_estimate ** 2) / snp_variance * (results.nobs / results.df_resid)

        # Adjust the standard errors
        std_raw = results.bse[v_name]
        std_adj = np.sqrt((std_raw ** 2) * ((self.total_obs - model_k) / (self.total_obs - (len(self.covariant) + 1))))

        # Use adjusted standard errors
        dist = stats.norm
        q = dist.ppf(1 - alpha / 2)
        lower_adj = estimate_adj - (q * std_adj)
        upper_adj = estimate_adj + (q * std_adj)

        # Calculate the two tailed chi squared test
        chi2tail = 1 - stats.chi2.cdf(estimate_adj, df=1)

        # Return the coefficient, standard errors, place values, obs, and lower + upper 95% CI
        return [estimate_adj, std_adj, results.pvalues[v_name], results.nobs, results.rsquared, chi2tail,
                lower_adj, upper_adj]
