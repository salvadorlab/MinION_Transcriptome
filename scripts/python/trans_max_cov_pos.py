import os
import sys
import pandas as pd

############################################
#
# identify the position with maximum number of covered reads for each transcript
# pandas required
# [trans_cov]: coverage per base for every transcript obtained using bedtools coverage
# [tran_max_cov]: max coverage of each transcript
#
# python trans_max_cov_pos.py [trans_cov] [tran_max_cov]
#
#############################################


trans_cov=sys.argv[1]
cov_pd = pd.read_csv(trans_cov, sep="\t", header=None)
trans_max_cov = cov_pd.groupby(3, as_index=False)[11].max() # group by transcript name, find max cov
trans_max_cov.iloc[:, 0] = trans_max_cov.iloc[:, 0].str[4:]
trans_max_cov.to_csv(sys.argv[2], sep="\t", index=False, header=False)


