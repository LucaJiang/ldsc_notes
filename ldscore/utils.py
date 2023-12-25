# Description: Utility functions for LD score regression
import numpy as np
import pandas as pd
import os, gzip, glob, logging
import multiprocessing as mp


#  Read z-score form .sumstats file
def read_sumstats(file_path):
    """
    Read z-score form .sumstats file
    :param file_path: path to .sumstats file
    :return: pandas dataframe with columns 'SNP' and 'Z'
    """
    df = pd.read_csv(file_path, sep="\t", header=0, usecols=["SNP", "Z"])
    logging.info("Read {} SNPs from {}".format(df.shape[0], file_path))
    ## delete rows without z-score
    df = df.dropna(subset=["Z"])
    logging.info("Remain {} SNPs with z-score".format(df.shape[0]))

    return df


# Merge sumstats with *.snplist file
# the remaining SNPs are those have r^2(L2)
def merge_sumstats_with_snplist(sumstats, snplist):
    """
    Merge sumstats with *.snplist file
    :param sumstats: pandas dataframe with columns 'SNP' and 'Z'
    :param snplist: pandas dataframe with columns 'SNP'
    :return: pandas dataframe with columns 'SNP' and 'Z'
    """
    merged = pd.merge(sumstats, snplist, on="SNP", how="inner")
    logging.info("Merged {} SNPs with *.snplist file".format(merged.shape[0]))

    return merged


# Get L2 of each SNP
## sub-function for parallel computing
def _merge_with_l2(ld_file, snplist):
    """
    Merge snplist with a single ldscore file
    :param ld_file: path to the ldscore file
    :param snplist: pandas dataframe with columns 'SNP' and 'Z'
    :return: merged pandas dataframe
    """
    with gzip.open(ld_file, "rt") as f:
        ldscore = pd.read_csv(f, sep="\t")
        return pd.merge(snplist, ldscore, on="SNP")


def get_l2(snplist, ld_dir):
    """
    Get L2 of each SNP
    :param snplist: pandas dataframe with columns 'SNP' and 'Z'
    :param ld_dir: path to LD score files
    :return: pandas dataframe with columns 'SNP', 'Z', 'CHR', 'BP', 'CM', 'L2'
    :algorithm: parallel {open LD score file -> merge with snplist}, then concat
    """
    ## get the dir of LD score files of all chromosomes
    ld_files = glob.glob(os.path.join(ld_dir, "*.l2.ldscore.gz"))

    with mp.Pool() as pool:
        merged_list = pool.starmap(
            _merge_with_l2, [(ld_file, snplist) for ld_file in ld_files]
        )
    logging.info("Merged {} LD score files".format(len(ld_files)))

    return pd.concat(merged_list)


# Calculate LD score of each SNP in a chromosome
## sub-function for parallel computing
def _calc_ldscore_by_chromosome(snplist, method="CM", window_size=1e-2):
    """
    Calculate LD score of each SNP in a chromosome
    :param snplist: pandas dataframe with columns 'SNP', '_distance_' and 'L2'
    :param method: method to calculate LD score, 'CM' or 'BP'
    :param window_size: window size to calculate LD score
    :return: pandas dataframe with columns 'SNP', 'CHR', 'Z', 'LDSCORE'
    :formula: $\ell = \sum_{i=1}^M \frac{l2_i}{M}$, where $M$ is the number of SNPs in the window
    :algorithm: for each snp, find the snps in the sliding window, then calculate the LD score
    """
    ## calculate LD score
    n = snplist.shape[0]
    ldscores = np.zeros(n)
    half_window_size = window_size / 2
    method_values = snplist[method].values
    l2_values = snplist["L2"].values
    start = 0
    end = -1
    window_sum = 0.0
    window_count = 0
    for i in range(n):
        while start < n and method_values[start] < method_values[i] - half_window_size:
            window_sum -= l2_values[start]
            window_count -= 1
            start += 1
        while (
            end + 1 < n
            and method_values[end + 1] <= method_values[i] + half_window_size
        ):
            end += 1
            window_sum += l2_values[end]
            window_count += 1
        ldscores[i] = window_sum / window_count
    snplist["LDSCORE"] = ldscores
    return snplist


def calc_ldscore(snplist, method="CM", window_size=1e-2):
    """
    Calculate LD score of each SNP, parallel by each chromosome
    :param snplist: pandas dataframe with columns 'SNP', 'CHR', '_distance_' and 'L2'
    :param method: method to calculate LD score, 'CM' or 'BP'
    :param window_size: window size to calculate LD score
    :return: pandas dataframe with columns 'SNP', 'Z', 'CHR', 'M', 'LDSCORE'; M is the number of SNPs in the window; where the LD score has been '/ M'
    """
    chromosomes = snplist["CHR"].unique()
    snplist["LDSCORE"] = np.nan
    with mp.Pool() as pool:
        ldscore_list = pool.starmap(
            _calc_ldscore_by_chromosome,
            [
                (snplist[snplist["CHR"] == chromosome], method, window_size)
                for chromosome in chromosomes
            ],
        )
    logging.info("Calculated LD score for {} chromosomes".format(len(chromosomes)))

    return pd.concat(ldscore_list)


# From sumstats file to ldscore file
def sumstats2ldsc(sumstats_file, reference_panel, N, method, window_size, out_file):
    """
    From sumstats file to ldscore file
    :param sumstats_file: path to sumstats file
    :param reference_panel: path to reference panel
    :param N: sample size of reference panel
    :param method: method to calculate LD score, 'CM' or 'BP'
    :param window_size: window size to calculate LD score
    :param out_file: path to output file
    :return: pandas dataframe with columns 'SNP', 'Z', 'CHR', 'M', 'LDSCORE'; M is the number of SNPs in the window, LDSCORE is N/M*\ell
    """
    ## get absolute path
    path_to_sumstats = os.path.abspath(sumstats_file)
    path_to_reference_panel = os.path.abspath(reference_panel)
    out_file = os.path.abspath(out_file)
    ## read sumstats
    sumstats = read_sumstats(path_to_sumstats)
    ## find *.snplist in reference panel
    path_to_snplist = glob.glob(os.path.join(path_to_reference_panel, "*.snplist"))[0]
    snplist = pd.read_csv(path_to_snplist, sep="\t", header=0, usecols=["SNP"])
    ## merge sumstats with snplist
    merged = merge_sumstats_with_snplist(sumstats, snplist)
    l2 = get_l2(merged, path_to_reference_panel)
    ## calculate LD score
    snp_ldscore = calc_ldscore(l2, method, window_size)
    snp_ldscore["LDSCORE"] *= N
    ## write to file
    snp_ldscore.to_csv(out_file, sep="\t", index=False)
    logging.info("Wrote LD score to {}".format(out_file))
    return snp_ldscore


#! TEST
if __name__ == "__main__":
    sumstats_file = "/Users/lucajiang/learn/CityU/ldsc_notes/data/full.sumstats"
    reference_panel = "/Users/lucajiang/learn/CityU/ldsc_notes/data/eur_w_ld_chr/"
    method = "CM"
    window_size = 1e-4
    logging.basicConfig(
        filename="test.log",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    ldscore = sumstats2ldsc(
        sumstats_file, reference_panel, 61220, method, window_size, "test.out"
    )
    print(ldscore.head())
