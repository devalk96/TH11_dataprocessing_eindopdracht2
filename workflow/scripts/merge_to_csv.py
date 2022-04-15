import sys
import pandas as pd
from scipy.io import arff


def main():
    # Get input
    infiles = snakemake.input
    outfile = snakemake.output.combined

    # Merge all files into dataframe
    df = pd.DataFrame()
    for file in infiles:
        ar = arff.loadarff(file)
        if df.empty:
            df = pd.DataFrame(ar[0])
        else:
            df = pd.concat([df, pd.DataFrame(ar[0])], axis=0)

    # Convert byte cols to str
    for col, dtype in df.dtypes.items():
        if dtype == object:
            df[col] = df[col].str.decode('utf-8')

    # Write to file
    df.to_csv(outfile, index=False)

    return 0


if __name__ == '__main__':
    sys.exit(main())
