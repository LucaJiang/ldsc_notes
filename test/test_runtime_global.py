import time, os, glob

start_time = time.time()

ld_dir = "data/eur_w_ld_chr/"
ld_files = glob.glob(os.path.join(ld_dir, "*.l2.ldscore.gz"))
ld_files = [ld_file for ld_file in ld_files if "old" not in ld_file]

end_time = time.time()

print("Execution time: {} seconds".format(end_time - start_time))
