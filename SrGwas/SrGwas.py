from miscSupports import validate_path, directory_iterator, load_yaml, FileOut, suppress_stdout
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

        # Set the gen file info
        self.gen_directory = self.args["path_to_gen_files"]
        self.gen_type = self.args["gen_type"]
        self.target_chromosome = self.args["target_chromosome"]

        # Set the output path for the memory files and load the file reference
        custom_meta_path(validate_path(self.args["memory_file_location"]))

        # Load the genetic reference, and sort both it and the external variables so they match on iid
        self.gen, self.variables,  = self._setup_variables()
        self.formula, self.phenotype = self._set_formula()

        # Set output file
        self.output = FileOut(validate_path(self.args["output_directory"]),
                              f"{self.args['output_name']}_Chr{self.target_chromosome}", "csv")
        self.output.write_from_list(["Snp"] + [iid for fid, iid in self.gen.iid])

        # Isolate which snps are to be used
        self.snp_ids = self._select_snps()

        # Start the method that has been assigned if method has been set
        if self.args["method"]:
            getattr(self, self.args["method"])()

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

        # Isolate the IID to match against the variables IID
        genetic_iid = [iid for fid, iid in gen.iid]

        sorting = []
        for index, iid in enumerate(variables["IID"].tolist()):
            # If the IID is in the genetic file, note its position in the genetic file for sorting
            if iid in genetic_iid:
                sorting.append(genetic_iid.index(iid))
            # Otherwise destroy this row
            else:
                variables.drop(index, axis=0, inplace=True)

        # Sort the dataframe so that the order of the dataframe is the same as the genetic file
        column_names = variables.columns
        variables = variables.to_numpy()
        variables = pd.DataFrame(variables[np.argsort(sorting)])
        variables.columns = column_names

        # Remove any IID in the genetic file that does not have variable information
        variable_iid = variables["IID"].tolist()
        gen = gen[[i for i, n in enumerate(genetic_iid) if n in variable_iid], :]

        for index, variable in enumerate(variables.columns):
            if index != self.args["variable_iid_index"]:
                variables[variable] = variables[variable].apply(pd.to_numeric)

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
        # todo Set id when id's are not random such as from summary stats
        raise NotImplementedError("Not yet in place")

    def create_genetic_residuals(self):

        for index, snp_i in enumerate(self.snp_ids):
            if index % 1000 == 0:
                print(f"{index} / {len(self.snp_ids)}")

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

            # Run the estimation, hiding its output as its unnecessary (akin to quietly)
            with suppress_stdout():
                results = ols_high_d_category(df, formula=f"Dosage~{self.formula}")

            # Extract the snp name and save the residuals
            snp_name = [self.gen.sid[snp_i].split(",")[1]]
            self.output.write_from_list(snp_name + results.resid.astype("string").tolist())
