from miscSupports import directory_iterator, chunk_list
from csvObject import CsvObject, write_csv
from pathlib import Path
import numpy as np


def main_call(out_dir, write_dir, headers):

    output = []
    for file in directory_iterator(out_dir):
        if ".log" not in file:
            csv_file = CsvObject(Path(output_dir, file))

            # Isolate the model values from the aggregated [snp] + [model 1, ... model N]
            for row in csv_file.row_data:
                snp, models = row[0], chunk_list(row[1:], len(headers))
                output.append([snp, models])

    print(f"For {len(output)} Snps")
    model_count = len(output[0][1])

    model_comp = []
    for i in range(model_count):
        print(f"For model {i+1}")

        # Write out the aggregated chromosome model data to a directory
        model_out = []
        for snp, model in output:
            model_out.append([snp] + model[i])
        write_csv(write_dir, f"Model{i + 1}", ["Snp"] + headers, model_out)

        # Append the comparision to a master list of models
        model_comp.append([f"Model {i+1}"] + [str(np.mean([float(values[vi]) for values in model_out]))
                                              for vi in range(1, 3)])

    # Write the model comp out
    write_csv(write_dir, "Model Comparision", ["Model", "Mean Coefficent", "Mean Standard Error", "Mean P Values"],
              model_comp)


if __name__ == '__main__':

    output_dir = r"C:\Users\Samuel\PycharmProjects\SR-GWAS\Testing\Output\OutputM5"
    write = r"Z:\UKB\GeographicID\Paper Data Extraction\SB_Papers\SW_GWAS\Output"
    model_headers = ["coef", "std_err", "pvalue", "obs", "95%lower", "95%upper"]

    main_call(output_dir, write, model_headers)
