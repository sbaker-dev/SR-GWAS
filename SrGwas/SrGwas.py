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

        # Set variable seek position on iid
        self.variables = validate_path(self.args["variables"])
        self.zipped = self.variables.suffix == ".gz"
        self.iid_location = self._set_seeks()

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

    def _set_seeks(self):
        """
        Exposure / external variables could be a vary large file, this allows us to index the expose file to the right
        location

        :return: A dict of {iid: seek}
        :rtype dict
        """

        iid_seeker = {}
        with open(self.variables, "rb") as file:

            # Skip header
            current_position = len(file.readline())

            # Note each iid location within the file
            for line_byte in file:
                iid = line_byte.decode().strip('\r\n').split(",")[self.args["variable_iid_index"]]
                iid_seeker[iid] = current_position
                current_position += len(line_byte)

            file.close()

        return iid_seeker

    def set_snp_ids(self):
        # todo Set id when id's are not random such as from summary stats
        raise NotImplementedError("Not yet in place")

    def create_genetic_residuals(self):
        print("Hello")
        print(self.gen.iid_count)
