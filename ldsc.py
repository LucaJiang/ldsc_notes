import argparse, logging, os, time, sys
from ldscore.calLDScore import sumstats2ldsc
import ldscore.irwls


# Input: sumstats_file, reference_panel, N, method, window_size, out_file
parser = argparse.ArgumentParser(
    description="Calculate heritability from summary statistics"
)
parser.add_argument(
    "--sumstats",
    "-s",
    type=str,
    required=True,
    help="Path to summary statistics file",
)
parser.add_argument(
    "--ref_panel",
    "-r",
    type=str,
    required=True,
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
    # Get args and Create Logger ------------------------------------------------
    args = parser.parse_args()
    start_time = time.time()
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Create a file handler
    file_handler = logging.FileHandler(os.path.join("results", args.out + ".log"))
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )

    # Create a stream handler
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )

    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    logging.info("Calculating LD score")
    variables = sumstats2ldsc(
        args.sumstats,
        args.ref_panel,
        args.N,
        args.method,
        args.window_size,
        args.out,
    )
    logging.info("Finished in {} seconds".format(time.time() - start_time))
    x = variables["LDSCORE"]
    y = variables["Z"] ** 2
