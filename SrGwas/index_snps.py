from miscSupports import validate_path, terminal_time, load_yaml
from csvObject import CsvObject, write_csv
from bgen_reader import custom_meta_path
from pysnptools.distreader import Bgen
from random import sample
from pathlib import Path


def set_snp_ids(yaml_path, gen_path, write_dir, file_name):
    """
    Isolate a subset of snps based on pre-defined named snps in a csv, passed as a str to snps_to_id, or a random
    set of snps of total == pre-defined int, where the int is set to snps_to_id.

    :param yaml_path: The path to the yaml file
    :type yaml_path: Path | str

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
    args = load_yaml(validate_path(yaml_path))
    custom_meta_path(validate_path(args["memory_file_location"]))
    gen = Bgen(str(validate_path(gen_path).absolute()))

    # If the snps_to_id is a str, assume a path and index the snps in the csv
    if isinstance(args["snps_to_id"], str):
        snps = CsvObject(validate_path(args["snps_to_id"]), set_columns=True)[0]
        snp_list = [f"{snp},{snp}" for snp in snps]
        snp_list = gen.sid_to_index(snp_list).tolist()

    # If an int, get a random sample of snps of population sid_count snp_to_id times.
    elif isinstance(args["snps_to_id"], int):
        snp_list = sample(range(gen.sid_count), args["snps_to_id"])

    # Else snp_to_id is invalid
    else:
        raise TypeError(f"set_snp_ids expects an int or a str yet was passed {type(args['snps_to_id'])}\n"
                        f"If you want to run a GWAS with all snps, simply set snps: null and ignore this method")

    write_csv(write_dir, f"{file_name}_Snps", ["Snp"], snp_list)
    print(f"Constructed snp id list {terminal_time()}")