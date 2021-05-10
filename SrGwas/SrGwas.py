from miscSupports import validate_path, directory_iterator, load_yaml, FileOut, suppress_stdout, terminal_time
from FixedEffectModel.api import ols_high_d_category
from bgen_reader import custom_meta_path
from pysnptools.distreader import Bgen
from csvObject import CsvObject
from pathlib import Path
import pandas as pd
import numpy as np
import re


class SrGwas:
    def __init__(self, args):

        # Load the args from the yaml file
        self.args = load_yaml(args)
        self.residual_run = self.args["residuals"]

        print(self.args)

        # Set the gen file info
        self.gen_directory = self.args["path_to_gen_files"]
        self.gen_type = self.args["gen_type"]
        self.target_chromosome = self.args["target_chromosome"]

        # Set the output path for the memory files and load the file reference
        custom_meta_path(validate_path(self.args["memory_file_location"]))
        self.logger = FileOut(self.args["output_directory"], f"{self.args['output_name']}_{self.target_chromosome}",
                              "log")

        self.logger.write(f"Setup {terminal_time()}")

        # Load the genetic reference, and sort both it and the external variables so they match on iid
        self.gen, self.variables = self._setup_variables()
        self.formula, self.phenotype = self._set_formula()
        self.logger.write(f"Set {self.gen.iid_count} in Genetic file and {len(self.variables)} in variable file.")

        # Set output file
        self.output = FileOut(validate_path(self.args["output_directory"]),
                              f"{self.args['output_name']}_Chr{self.target_chromosome}", "csv")
        if self.residual_run:
            self.output.write_from_list(["Snp"] + [iid for fid, iid in self.gen.iid])
        else:
            self.output.write_from_list(["Snp"] + ["coef", "std_err", "pvalue", "95%lower", "95%upper"])

        # Start the method that has been assigned if method has been set
        if self.args["method"]:
            getattr(self, self.args["method"])()

    def __repr__(self):
        return f"SrGwas object Controller"

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

    def _setup_variables(self):
        """
        The order of IID in genetic file may not equal to submission, this sorts the arrays to be equivalent.

        :return: Bgenfile for this chromosome as well as a pandas dataframe of the external variables
        :rtype: (Bgen, pd.DataFrame)
        """

        # Load the variables as pandas dataframe and setup the reference genetic file for this chromosome
        variables = pd.read_csv(validate_path(self.args["variables"]))
        gen = Bgen(self._select_file_on_chromosome())
        self.logger.write(f"...Loaded external variables {terminal_time()}")

        # Recast IID as an int
        variables["IID"] = [int(re.sub(r'[\D]', "", iid)) for iid in variables["IID"].tolist()]

        # Isolate the IID to match against the variables IID and create the reference
        genetic_iid = np.array([int(re.sub(r'[\D]', "", iid)) for fid, iid in gen.iid])
        genetic_position = gen.iid

        # Remove any IID that is in the external data array but not in the genetic array
        out = np.in1d(variables["IID"].to_numpy(), genetic_iid)
        variables = variables[out]

        # Remove any IID that is in the genetic array but not in the external data
        out = np.in1d(genetic_iid, variables["IID"].to_numpy())
        genetic_iid = genetic_iid[out]
        genetic_position = genetic_position[out]

        # Sort both arrays to be in the same order
        variables = variables.sort_values(by=['IID'], ascending=True)
        gen = gen[gen.iid_to_index(genetic_position[np.argsort(genetic_iid)]), :]

        for index, variable in enumerate(variables.columns):
            if index != self.args["variable_iid_index"]:
                variables[variable] = variables[variable].apply(pd.to_numeric)

        self.logger.write(f"Setup external reference {terminal_time()}")
        return gen, variables

    def _set_formula(self):
        """
        Set the stats model left hand side formula for this run

        :return: The string formula for the left hand side of the equation
        :rtype: (str, str)
        """
        # Validate each variable type
        [self._validate_variable(cont, "Continuous") for cont in self.args["continuous_variables"]]
        [self._validate_variable(cont, "Fixed_Effect") for cont in self.args["fixed_effect_variables"]]
        [self._validate_variable(cont, "Cluster") for cont in self.args["cluster_variables"]]

        # Validate the phenotype if set (won't necessarily need it for residual runs)
        if self.args["phenotype"]:
            self._validate_variable(self.args["phenotype"], "Phenotype")

        # Join Each variable with a +
        cont = "+".join([cont for cont in self.args["continuous_variables"]])
        fe = "+".join([cont for cont in self.args["fixed_effect_variables"]])
        cl = "+".join([cont for cont in self.args["cluster_variables"]])

        if len(cont) == 0:
            raise IndexError("No Variables provided to continuous_variables")

        if len(self.args["fixed_effect_variables"]) == 0 and len(self.args["cluster_variables"]) == 0:
            return cont, self.args["phenotype"]

        elif len(self.args["fixed_effect_variables"]) > 0 and len(self.args["cluster_variables"]) == 0:
            return f"{cont}|{fe}", self.args["phenotype"]

        elif len(self.args["fixed_effect_variables"]) > 0 and len(self.args["cluster_variables"]) > 0:
            return f"{cont}|{fe}|{cl}", self.args["phenotype"]

        elif len(self.args["fixed_effect_variables"]) == 0 and len(self.args["cluster_variables"]) > 0:
            raise IndexError("Invalid formula, clustering can only be undertaken when FE specified.")

        else:
            raise IndexError(f"Unknown formula specification of lengths {len(cont)}, {len(fe)}, {len(cl)}")

    def _validate_variable(self, v, var_type):
        """Check the variable exists within the columns"""
        assert v in self.variables.columns, f"{var_type} variable {v} not in variables {self.variables.columns}"

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
        print(self.gen.sid_count)

    def gwas(self):
        """
        Create genetic residuals by regressing your covariant on the snp or run a more traditional gwas of

        phenotype ~ dosage + covariant_1 + ... covariant_N

        :return: Nothing, write line to fine when residuals have been estimated
        :rtype: None
        """
        # Isolate which snps are to be used
        snp_ids = self._select_snps()

        for index, snp_i in enumerate(snp_ids):
            if index % 1000 == 0:
                print(f"{index} / {len(snp_ids)}")
                self.logger.write(f"{index} / {len(snp_ids)}")

            # Instance the memory for all individuals (:) for snp i
            current_snp = self.gen[:, snp_i]

            # Transform bgen dosage of [0, 1, 0] -> 0, 1, or 2 respectively.
            dosage = sum(np.array([snp * i for i, snp in enumerate(current_snp.read(dtype=np.int8).val.T)],
                                  dtype=np.int8))[0]

            # Combine the dataset for regression
            dosage = pd.DataFrame(dosage)
            dosage.columns = ["Dosage"]
            data = [pd.DataFrame(dosage), self.variables]
            df = pd.concat(data, axis=1)

            snp_name = [self.gen.sid[snp_i].split(",")[1]]
            if self.residual_run:
                # Run the estimation, hiding its output as its unnecessary (akin to quietly)
                with suppress_stdout():
                    results = ols_high_d_category(df, formula=f"Dosage~{self.formula}")

                # Extract the snp name and save the residuals
                self.output.write_from_list(snp_name + results.resid.astype("string").tolist())

            else:
                with suppress_stdout():
                    results = ols_high_d_category(df, formula=f"{self.phenotype}~Dosage+{self.formula}")

                # Extract the regression results
                regression_results = [
                    results.params["Dosage"],
                    results.bse["Dosage"],
                    results.pvalues["Dosage"]] + results.conf_int().loc["Dosage"].tolist()
                regression_results = [str(r) for r in regression_results]

                self.output.write_from_list(snp_name + regression_results)
