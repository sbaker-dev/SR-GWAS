from miscSupports import z_scores, flip_list
from random import uniform, randint
from csvObject import write_csv


if __name__ == '__main__':

    # Setup
    sample_size = 483
    pc_count = 10
    example_snp = 10
    output_directory = r"C:\Users\Samuel\PycharmProjects\SR-GWAS\Data"

    # Setup IID
    iid = [f"sample_{i}" for i in range(sample_size)]

    # Setup basic Identifiers
    gender = [randint(0, 1) for _ in range(sample_size)]
    yob = [randint(1934, 1971) for _ in range(sample_size)]

    # Phenotype of BMI
    bmi = [uniform(14.4, 30.2) for _ in range(sample_size)]
    output = [iid, bmi, gender, yob]

    # Add PC's then write as covariant file
    for i in range(pc_count):
        output.append(z_scores([randint(0, 1000) for _ in range(sample_size)]))
    headers = ["IID", "BMI", "Gender", "Age"] + [f"PC{i}" for i in range(1, 11)]
    write_csv(output_directory, "Covariant", headers, flip_list(output))

    # Add Example Snps then write as covariant + snp file
    for i in range(example_snp):
        output.append([randint(0, 2) for _ in range(sample_size)])
    headers = headers + [f"rs{i}{i+1}{i+2}" for i in range(example_snp)]
    write_csv(output_directory, "CovariantSnp", headers, flip_list(output))
