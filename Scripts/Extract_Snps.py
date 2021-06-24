from miscSupports import open_setter, decode_line
from random import sample, seed
from operator import itemgetter
from csvObject import write_csv
from pathlib import Path


if __name__ == '__main__':

    p_values = []

    root = Path(r"Z:\UKB\GeographicID\Paper Data Extraction\SB_Papers\SW_GWAS\SNP_gwas_mc_merge_nogc.tbl (1).uniq.gz")
    with open_setter(root)(root) as f:

        for index, line in enumerate(f):
            if index % 100000 == 0:
                print(f"Extracted {index} snps")

            line_decoded = decode_line(line, True, "\t")

            p_values.append([line_decoded[-2], line_decoded[0]])

    # Set random seed for replication
    seed(10)
    random_snps = sample(range(len(p_values)), 10000)
    largest_hits = [[snp] for p, snp in sorted(p_values, key=itemgetter(0))[:10000]]
    random_snps = [[p_values[i][1]] for i in random_snps]

    # Write snp files
    write_csv(r"Z:\UKB\GeographicID\Paper Data Extraction\SB_Papers\SW_GWAS", "SNP_Random", ["Snp"], random_snps)
    write_csv(r"Z:\UKB\GeographicID\Paper Data Extraction\SB_Papers\SW_GWAS", "SNP_Top", ["Snp"], random_snps)