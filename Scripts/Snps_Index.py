from SrGwas.index_snps import set_snp_ids
from miscSupports import terminal_time
from pathlib import Path

if __name__ == '__main__':

    gen_dir = r"/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen"
    mem_loc = r"/mnt/storage/scratch/ow18390/pyGenic/meta_data"

    random = r"/mnt/storage/scratch/ow18390/SW_GWAS/SNP_Random.csv"
    top_10k = r"/mnt/storage/scratch/ow18390/SW_GWAS/SNP_Top.csv"

    out_random = r"/mnt/storage/scratch/ow18390/SW_GWAS/Random_Snps"
    out_top = r"/mnt/storage/scratch/ow18390/SW_GWAS/Top_Snps"

    for i in range(1, 23):

        gen_path = Path(gen_dir, f"data.chr{str(i).zfill(2)}.bgen")
        print(f"Loading {gen_path}: {terminal_time()}")

        set_snp_ids(mem_loc, random, gen_path, out_random, f"SnpsChr{i}")
        set_snp_ids(mem_loc, top_10k, gen_path, out_top, f"SnpsChr{i}")
