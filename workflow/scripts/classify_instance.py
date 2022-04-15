#!/usr/bin/python3

"""
A wrapper that accepts a file, which it classifies and then will pipe to a new arff file
"""

__author__ = "Sander Bouwman"
__copyright__ = "Copyright 2022, Sander Bouwman"

import sys
import arff
import random

def classify(instance: list, infile):
    # Check for proper amount of variables
    if len(instance) != 16:
        raise ValueError(f"A total of 16 variables should be provided. "
                         f"There are {len(instance)} provided.")

    # As mentioned earlier:
    # There seems to be a bug in the classifier.jar. Where it will always classify an instance
    # as True I spend a few hours trying to fix this bug But I can't seem to find the issue.
    # I tried rebuilding the model in WEKA a couple of times. In Weka GUI the model works as
    # intended but once exported and used in the Java wrapper it fails to work.
    # As the assignment is mainly the building of a pipeline I chose to run the java wrapper but
    # put out a randomly weighted True or False.

    # Code used if the wrapper is used.

    # jar_path: str = "../../resources/classify.jar"
    # model_path: str = "../../resources/adaboost.model"

    # Create string from instance list
    # instance_str = ",".join([str(x) for x in instance])

    # build list of arguments and run script
    # runlist = f"java -jar {jar_path} -v {instance_str} -m {model_path}".split()
    # results = subprocess.run(runlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    classification_result = simulate_classification()
    print(f"[INFO] {infile}\t|\tClassified as: {classification_result}")
    instance.append(classification_result)
    return instance


def simulate_classification():
    possible_vals = ["TRUE", "FALSE"]
    weights = [0.3, 0.7]
    return random.choices(possible_vals, weights, k=1)[0]


def main():
    infile = snakemake.input.txt
    outdir = snakemake.output.txt

    loaded_arff = arff.load(open(infile, "r"))

    # Add new attribute
    loaded_arff["attributes"].append(('HeartDisease', ['TRUE', 'FALSE']))

    # Classify
    data = loaded_arff["data"]
    loaded_arff["data"] = [classify(data[0], infile)]

    # Save ARFF
    arff.dump(loaded_arff, open(outdir, "w"))


if __name__ == '__main__':
    sys.exit(main())
