#!/usr/bin/python3

"""
Splits out a received arff file into arff files containing only a single instance.
"""
__author__ = "Sander Bouwman"
__copyright__ = "Copyright 2022, Sander Bouwman"

import os
import sys
import copy
import arff


def main():
    infile = snakemake.input[0]
    outdir = snakemake.output[0]

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    arr = arff.load(open(infile, "r"))

    for i, _ in enumerate(arr["data"]):
        new_arr = copy.deepcopy(arr)
        new_arr["data"].clear()
        new_arr["data"].append(arr["data"][i])
        arff.dump(new_arr, open(f"{outdir}/{i}.arff", "w"))


if __name__ == '__main__':
    sys.exit(main())
