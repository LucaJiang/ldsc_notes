import cProfile
from ldscore.calLDScore import sumstats2ldsc
import pstats, sys


def main():
    sumstats_file = "/Users/lucajiang/learn/CityU/ldsc_notes/data/full.sumstats"
    reference_panel = "/Users/lucajiang/learn/CityU/ldsc_notes/data/eur_w_ld_chr/"
    method = "CM"
    window_size = 1e-4

    ldscore = sumstats2ldsc(
        sumstats_file,
        reference_panel,
        method,
        window_size,
        "results/test",
    )


if __name__ == "__main__":
    cProfile.run("main()", "test/output.prof")

    # dump stats to a .txt file
    with open("test/output.txt", "w") as f:
        p = pstats.Stats("test/output.prof", stream=f)
        p.sort_stats("cumulative")

        # print top 10 functions by cumulative time
        p.print_stats(10)
        f.write("\n")

        # print callers of top 10 functions
        p.print_callers(10)
        f.write("\n")

        # print callees of top 10 functions
        p.print_callees(10)
