import argparse, logging, os, time, sys
from ldscore.calLDScore import sumstats2ldsc
import ldscore.irwls as irwls
import pandas as pd
import numpy as np


# Input: sumstats_file, reference_panel, N, method, window_size, out_file
parser = argparse.ArgumentParser(
    description="Calculate LD score or heritability from summary statistics"
)
parser.add_argument(
    "--mission",
    "-M",
    type=str,
    default="all",
    help="Mission to run, 'all', 'ldsc', or 'h2'",
)
parser.add_argument(
    "--sumstats",
    "-s",
    type=str,
    help="Path to summary statistics file",
)
parser.add_argument(
    "--ref_panel",
    "-r",
    type=str,
    help="Path to reference panel",
)
parser.add_argument(
    "--N",
    "-n",
    type=int,
    default=61220,
    help="Sample size of reference panel",
)
parser.add_argument(
    "--method",
    "-m",
    type=str,
    default="CM",
    help="Method to calculate LD score, 'CM' or 'BP'",
)
parser.add_argument(
    "--window_size",
    "-w",
    type=float,
    default=1e-2,
    help="Window size to calculate LD score",
)
parser.add_argument(
    "--out",
    "-o",
    type=str,
    default="run",
    help="Name of output file",
)

if __name__ == "__main__":
    start_time = time.time()

    # Get args and Create Logger ------------------------------------------------
    args = parser.parse_args()
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Create a file handler
    file_handler = logging.FileHandler(args.out + ".log")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S"
        )
    )

    # Create a stream handler
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S"
        )
    )

    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    logging.info("Mission: {}".format(args.mission))
    logging.info("Command: {}".format(" ".join(sys.argv)))

    # Run Model ---------------------------
    if args.mission == "ldsc" or args.mission == "all":
        logging.info("Calculating LD score")
        variables = sumstats2ldsc(
            args.sumstats,
            args.ref_panel,
            args.method,
            args.window_size,
            args.out,
        )
        logging.info("LD score: Finished in %.2f seconds" % (time.time() - start_time))
    if args.mission == "h2":
        # only calculate heritability, need to read parameters from file
        variables = pd.read_csv(args.sumstats, sep="\t")

    if args.mission == "all" or args.mission == "h2":
        logging.info("Calculating heritability")
        M = variables.shape[0]
        x = variables["LDSCORE"].values.reshape(-1, 1) * args.N / M
        y = variables["Z"].values.reshape(-1, 1) ** 2

        #! TEST
        coef_df = pd.DataFrame({"LDSCORE": x.flatten(), "Z^2": y.flatten()})
        coef_df.to_csv(args.out + "_coef.txt", sep="\t", index=False)
        logging.info("Wrote coefficients to {}".format(args.out + "_coef.txt"))

        x = np.concatenate([np.ones_like(x), x], axis=1)  # add a column of intercept
        print(x.shape, y.shape)  #! TEST
        irwls = irwls.IRLS(x, y)
        irwls.regression()
        reg_intercept = irwls.get_intercept()
        reg_coefficients = irwls.get_coefficients()
        logging.info("Intercept: %.4f" % reg_intercept)
        logging.info("h^2: %.4f" % reg_coefficients)
        logging.info("H^2: Finished in {} seconds".format(time.time() - start_time))

    logging.info("Mission completed.\n\n")
