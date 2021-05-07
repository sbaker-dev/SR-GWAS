from FixedEffectModel.api import ols_high_d_category
from pysnptools.distreader import Bgen
from miscSupports import validate_path, terminal_time
from contextlib import contextmanager
from random import randint
from pathlib import Path
import pandas as pd
import numpy as np
import sys
import os


@contextmanager
def suppress_stdout():
    """
    Suppress package print statements

    Note
    ----
    FixedEffectModels prints every step, for every regression. Might desirable on a single regression where the user
    will read the terminal but in our case it just clutters .log's with millions of lines of output

    Source: https://stackoverflow.com/questions/2125702/how-to-suppress-console-output-in-python
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


class ResOut:
    def __init__(self, write_directory, write_name):
        """
        Working with 100k's of results will eat memory if we store residual until final complication. As such we want to
        write the csv's files dynamically

        :param write_directory: Directory of output file
        :type write_directory: Path | str

        :param write_name: File name
        :type write_name: str
        """

        self.file = open(Path(write_directory, f"{write_name}.csv"), "w")

    def write(self, line):
        """Write and flush a line to a file"""
        self.file.write(f"{line}\n")
        self.file.flush()

    def write_from_list(self, values_list):
        """Write a list of variables to a file as a comma delimited line"""
        self.file.write(f"{','.join(values_list)}\n")

    def close(self):
        """Close the log via objected rather than attribute"""
        self.file.close()


if __name__ == '__main__':

    # Set the path to the example file you downloaded here as well as the output directory. Do not remove the forward
    # 'r' as this allows the string to be read as literal
    path_to_gen_file = r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\ByChromosome\EUR.ldpred_1.bgen"
    output_directory = r"I:\Work\Genetics\Residuals"

    # Load gen file
    # It will take longer the first time you use it as it has to create pysnptools metadat .mmm files, which act as a
    # faster .bim /.bgi equivalent. Finalised system will also utilise pyGenicParser which is an extended version of
    # pyBgen written by myself that can use .bim in lue of .mmm files when hard-drive / scratch space is at a premium.
    # Will also be able to use .bed/bim/fam files for plink interfacing. See snipped of the pyGenicPipeline i am the
    # author of at line 212 method isolate_raw_snps of the following for reference if required.
    # https://github.com/sbaker-dev/pyGenicPipeline/blob/main/pyGenicPipeline/core/Loaders/commonGenetic.py
    gen = Bgen(validate_path(path_to_gen_file))

    # Create some dummy values for gender and district Fixed effects, iid count of example dataset is 483 so district
    # FE set lower than actual so that there is within district groups.
    gender = [randint(0, 1) for _ in range(gen.iid_count)]
    district = [randint(0, 50) for _ in range(gen.iid_count)]
    print(f"Total number of individuals {gen.iid_count}")

    # Output file takes for form of a matrix of M x N in csv form, where M rows represent the number of snps that have
    # residuals for N number of individuals. As a note, csv files are unlikely to be sufficient for storing this amount
    # of data, but we can use the plink / bgen IoStream byte code write/unpack logic if we generally need this amount
    # of data
    file = ResOut(validate_path(output_directory), "ResidualsOVERIDE")

    # Construct IO stream to write out to and write the header of Snp + [IID1, IID2, ... IID(N)].
    # Bgen files store [variant id, rs_id], we just want the rs_id hence the [1]; see https://bit.ly/2J0C1kC
    # Loader will utilise load_variants of the same file as referenced in load gen line 158 which makes it bgen, bed
    # and non .mmm file compliant.
    file.write_from_list(["Snp"] + [iid for fid, iid in gen.iid])

    # For each snp in the current chromosome file
    for i in range(gen.sid_count):

        if i % 100 == 0:
            print(f"{i}/{gen.sid_count}: {terminal_time()}")

        # Instance the memory for all individuals (:) for snp i
        current_snp = gen[:, gen.sid_to_index([gen.sid[i]])]

        # Transform bgen dosage of [0, 1, 0] -> 0, 1, or 2 respectively.
        dosage = sum(np.array([snp * i for i, snp in enumerate(current_snp.read(dtype=np.int8).val.T)], dtype=np.int8))

        # Set the dosage data into an array with covariant's
        df = pd.DataFrame([dosage[0]]).T
        df.columns = ["Dosage"]
        df["Gender"] = gender
        df["District"] = district

        # Formula takes for the form of 'dependent variable ~ continuous variable|fixed_effect|clusters'
        # In this case we are not clustering but can do. Variables that are considered as fixed_effects are de-meaned.
        # Python package of https://pypi.org/project/FixedEffectModel/ aims to replicated that stata equivalent of
        # high dimension fixed effects transformations/modeling of reghdfe of http://scorreia.com/software/reghdfe/
        # which i use myself in stata
        formula = "Dosage~Gender|District"

        with suppress_stdout():
            # Run the estimation, hiding its output as its unnecessary (akin to quietly)
            results = ols_high_d_category(df, formula=formula)

        # Bgen files store snps as [variant id, rs_id], we just want the rs_id hence the [1]; see https://bit.ly/2J0C1kC
        snp_name = [gen.sid[i].split(",")[1]]
        file.write_from_list(snp_name + results.resid.astype("string").tolist())

    print(f"Finished time {terminal_time()}")
