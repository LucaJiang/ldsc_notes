# Calculate LD score from sumstats file
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
    df = pd.read_csv(file_path, sep="\t", header=0, usecols=["SNP", "Z", "A1", "A2"])
    logging.info("Read {} SNPs from {}".format(df.shape[0], file_path))
    ## delete rows without z-score
    df.dropna(subset=["Z"], inplace=True)
    # check if there are duplicated SNPs
    duplicated = df.duplicated(subset=["SNP"])
    if np.any(duplicated):
        df = df[~duplicated]
        logging.info("Removed {} duplicated SNPs".format(np.sum(duplicated)))
    logging.info("Remain {} unique SNPs with z-score".format(df.shape[0]))

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
    # check A1 and A2 match
    # order is not important
    sumstats[["A1", "A2"]] = np.sort(sumstats[["A1", "A2"]].values, axis=1)
    snplist[["A1", "A2"]] = np.sort(snplist[["A1", "A2"]].values, axis=1)

    merged = pd.merge(sumstats, snplist, on=["SNP", "A1", "A2"], how="inner")
    logging.info(
        "Remain {} SNPs after merged with *.snplist file".format(merged.shape[0])
    )
    return merged[["SNP", "Z"]]  # only keep 'SNP' and 'Z'


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
    ldscore.set_index("SNP", inplace=True)
    return ldscore.join(snplist, on="SNP", how="inner")


def get_l2(snplist, ld_dir):
    """
    Get L2 of each SNP
    :param snplist: pandas dataframe with columns 'SNP' and 'Z'
    :param ld_dir: path to LD score files
    :return: pandas dataframe with columns 'SNP', 'Z', 'CHR', 'BP', 'CM', 'L2'
    :algorithm: parallel {open LD score file -> merge with snplist}, then concat
    """
    chromosomes_id = np.arange(1, 23)
    ## get the dir of LD score files of all chromosomes
    ld_files = [
        os.path.join(ld_dir, "{}.l2.ldscore.gz".format(i)) for i in chromosomes_id
    ]
    logging.info("Merged {} LD score files".format(len(ld_files)))

    ## set SNP as index to accelerate
    snplist.set_index("SNP", inplace=True)
    if False:  # MPI
        with mp.Pool() as pool:
            merged_list = pool.starmap(
                _merge_with_l2, [(ld_file, snplist) for ld_file in ld_files]
            )
        # concatenate all chromosomes
        merged_list = pd.concat(merged_list)
    else:
        merged_list = pd.concat(
            [_merge_with_l2(ld_file, snplist) for ld_file in ld_files]
        )
    len_before = merged_list.shape[0]
    # delete MAF <= 0.01
    merged_list = merged_list[merged_list["MAF"] > 0.01]
    logging.info(
        "Removed {} SNPs with MAF <= 0.01".format(len_before - merged_list.shape[0])
    )

    return merged_list


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
    # estimate the end of the window to accelerate
    end = half_window_size / (method_values[1] - method_values[0] + 1e-8) * 0.8
    end = int(end)
    window_sum = np.sum(l2_values[: end + 1])
    for i in range(n):
        while start <= i and method_values[start] < method_values[i] - half_window_size:
            window_sum -= l2_values[start]
            start += 1
        while (
            end + 1 < n
            and method_values[end + 1] <= method_values[i] + half_window_size
        ):
            end += 1
            window_sum += l2_values[end]
        ldscores[i] = window_sum / (end - start + 1)
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
    chromosomes.sort()
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
def sumstats2ldsc(**kwargs):
    """
    From sumstats file to ldscore file
    :return: pandas dataframe with columns 'SNP', 'Z', 'CHR', 'M', 'l2'; M is the number of SNPs in the window, LDSCORE is l2
    """
    ## get absolute path
    path_to_sumstats = os.path.abspath(kwargs["sumstats"])
    path_to_reference_panel = os.path.abspath(kwargs["ref_panel"])
    out_file = os.path.abspath(kwargs["out"])
    ## read sumstats
    sumstats = read_sumstats(path_to_sumstats)
    ## find *.snplist in reference panel
    path_to_snplist = glob.glob(os.path.join(path_to_reference_panel, "*.snplist"))[0]
    snplist = pd.read_csv(path_to_snplist, sep="\t", header=0)
    ## merge sumstats with snplist
    merged = merge_sumstats_with_snplist(sumstats, snplist)
    l2 = get_l2(merged, path_to_reference_panel)
    if kwargs["method"] is not None and kwargs["window_size"] is not None:
        ## calculate LD score within window measure by method
        snp_ldscore = calc_ldscore(l2, kwargs["method"], kwargs["window_size"])
    else:
        snp_ldscore = l2
    ## write to file
    snp_ldscore.to_csv(out_file + ".txt", sep="\t", index=False)
    logging.info("Wrote LD score to {}".format(out_file))
    return snp_ldscore


#! TEST
import sys

if __name__ == "__main__":
    sumstats_file = "./data/full.sumstats"
    reference_panel = "./data/eur_w_ld_chr/"
    method = "CM"
    window_size = 1e-4

    dict = {
        "sumstats": sumstats_file,
        "ref_panel": reference_panel,
        "method": method,
        "window_size": window_size,
        "out": "results/test",
    }

    ldscore = sumstats2ldsc(**dict)
    print(ldscore.head())
