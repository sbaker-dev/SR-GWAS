from miscSupports import validate_path, terminal_time
from csvObject import CsvObject, write_csv
from bgen_reader import custom_meta_path
from pysnptools.distreader import Bgen
from pathlib import Path


def set_snp_ids(memory_location, snps_to_id, gen_path, write_dir, file_name):
    """
    Isolate a subset of snps based on pre-defined named snps in a csv, passed as a str to snps_to_id, or a random
    set of snps of total == pre-defined int, where the int is set to snps_to_id.

    :param memory_location: Location of bgen memory file
    :type memory_location: Path | str

    :param snps_to_id: Location of snps csv to id
    :type snps_to_id: Path | str

    :param gen_path: The path to the genetic file
    :type gen_path: Path | str

    :param write_dir: The directory to write the snp index csv file to
    :type write_dir: Path | str

    :param file_name: The name of the snp index file
    :type file_name: str

    :return: Nothing, write the id's to a csv then stop
    :rtype: None

    :raise TypeError: If a str / int is not passed
    """

    # Load the args dict, then set the custom write location for the bgen file memory files and load the genetic ref
    custom_meta_path(validate_path(memory_location))
    gen = Bgen(str(validate_path(gen_path).absolute()))

    # If the snps_to_id is a str, assume a path and index the snps in the csv
    snps = CsvObject(validate_path(snps_to_id), set_columns=True)[0]
    snp_list = [f"{snp},{snp}" for snp in snps]
    snp_list = gen.sid_to_index(snp_list).tolist()

    # Write the snps out
    write_csv(write_dir, f"{file_name}_Snps", ["Snp"], snp_list)
    print(f"Constructed snp id list {terminal_time()}")
