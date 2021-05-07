from miscSupports import validate_path, directory_iterator, load_yaml
from bgen_reader import custom_meta_path
from pysnptools.distreader import Bgen
from pathlib import Path
import re


class SrGwas:
    def __init__(self, args):

        self.args = load_yaml(args)

        # Set the gen file info
        self.gen_directory = self.args["path_to_gen_files"]
        self.gen_type = self.args["gen_type"]
        self.target_chromosome = self.args["target_chromosome"]

        # Set the output path for the memory files
        custom_meta_path(validate_path(self.args["memory_file_location"]))

        self.gen_file = Bgen(self._select_file_on_chromosome())

        print(self.gen_file.iid_count)

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
