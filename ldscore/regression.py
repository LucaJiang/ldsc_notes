n_block = 200

import numpy as np
import pandas as pd

coef = pd.read_csv("./results/test_coef.txt", sep="\t")
x = coef["LDSCORE"].values.reshape(-1, 1)
y = coef["Z^2"].values.reshape(-1, 1)
chisq = y
