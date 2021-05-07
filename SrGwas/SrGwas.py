from miscSupports import validate_path, directory_iterator, load_yaml, FileOut
from bgen_reader import custom_meta_path
from pysnptools.distreader import Bgen
from csvObject import CsvObject
from pathlib import Path
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
        self.gen = Bgen(self._select_file_on_chromosome())

        # Set output file
        output = FileOut(validate_path(self.args["output_directory"]), self.args["output_name"], "csv")
        output.write_from_list(["Snp"] + [iid for fid, iid in self.gen.iid])

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

    def _select_snps(self):
        """
        We may only want to run a subset of snps. If so, then this loads the snp indexes from a csv. Else, just return
        all the snp ids

        :return: A list of snp ids
        :rtype: list[snp]
        """
        if self.args["snps"]:
            return CsvObject(validate_path(self.args["snps"]), set_columns=True)[0]
        else:
            return [i for i in range(self.gen.sid_count)]

    def set_snp_ids(self):
        # todo Set id when id's are not random such as from summary stats
        raise NotImplementedError("Not yet in place")

    def create_genetic_residuals(self):
        print("Hello")
