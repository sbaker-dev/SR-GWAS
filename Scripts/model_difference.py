from miscSupports import directory_iterator, flatten
from collections import Counter
from csvObject import CsvObject
from pathlib import Path
import numpy as np
import re


def main_call(out_dir):

    all_coefficients = []
    all_se = []
    for i in range(1, 23):
        model_csv = []
        for model_dir in directory_iterator(out_dir, False):
            model_csv.append(CsvObject(isolate_chr_file(out_dir, model_dir, i), set_columns=True))

        unique_snps = Counter(flatten([csv[0] for csv in model_csv]))
        snp = set([snp for snp, count in unique_snps.items() if count == 4])
        model_rows = [[row for row in model.row_data if row[0] in snp] for model in model_csv]

        all_coefficients.append([[float(row[1]) for row in model] for model in model_rows])
        all_se.append([[float(row[2]) for row in model] for model in model_rows])

    model_info = []
    for i in range(4):
        coefficients = flatten([chromosome[i] for chromosome in all_coefficients])
        se = flatten([s[i] for s in all_se])
        model_info.append([np.mean(coefficients), np.mean(se)])

    for r in model_info:
        print(r)


def isolate_chr_file(out_dir, model_dir, chromosome):
    for file in directory_iterator(Path(out_dir, model_dir)):
        if ".log" not in file:
            if chromosome == int(re.sub(r'[\D]', "", file.split("_")[1])):
                return Path(out_dir, model_dir, file)


if __name__ == '__main__':

    output_dir = r"C:\Users\Samuel\PycharmProjects\SR-GWAS\Testing\Output"

    main_call(output_dir)
