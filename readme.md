# Fast classification using a Snakemake pipeline

## Table of content
- [Introduction](#Introduction)
- [Pipeline](#Pipeline)
    * [Requirements](#Requirements)
      * [Python](#Python-3.5+)
      * [R](#R)
    * [Pipeline steps](#Pipeline-steps)
    * [Running the pipeline](#Running-the-pipeline)
- [Directory structure](#Directory-structure)
- [Contact](#contact)

## Introduction
Currently, heart disease is the number one leading cause of death, with more than 17 million 
deaths annually. These deaths account for more than 30% of all deaths worldwide. 
An early detection of possible heart disease could prevent a lot of future deaths. 
Machine learning could play a significant role in the early detection of heart disease. 
Currently, heart disease is diagnosed by a combination of blood test, ECG, breathing tests 
and chest x-rays test. For mass screening, it would be best to provide an accurate prediction 
of developing heart disease with the least amount of tests. Preferable, 
tests that can be performed at a general practitioner, which prevents unnecessary hospital visits. 
The goal of this pipeline was to create an accurate model that uses easy data which could 
easily be gathered at a general practitioner.

## Pipeline 
This pipeline is created for the purpose of predicting heart disease in patients.
The input file of the pipeline is an ARFF file. This file is not completely clean as is contains 
instances where there are missing values. The classifier can not handle missing values so these 
instances will be removed from the dataset.

### Requirements
#### Python 3.5+
| Name      | Version |
|-----------|---------|
| snakemake | 7.3.8   |
| pandas    | 1.4.2   |
| pandas    | 1.4.2   |
| scipy     | 1.8.0   |
| arff      | 0.9     |
  
#### R
| Name    | Version |
|---------|---------|
| ggplot2 | 3.3.5   |


### Pipeline steps

![Pipeline](https://github.com/devalk96/T11_Dataprocessing_Eindopdracht/blob/master/workflow/dag.png)

#### Step 1
##### Clean *ARFF* file from missing values
Uses script: *clean_arff.R*  
The pipeline uses arff files as input. Rows with missing values will be removed as model doesn't 
accept missing values

#### Step 2
##### Split arff files out with *split_arff.py*
For parallelization purposes it is better to split out the arff files into arff files which only
contains a single instance.

#### Step 3
##### Classify using *classify_instance.py*
Classify instance calls *classifier.jar* (located in resources.) This model will classify for 
heart disease.

#### Step 4
##### Merge all the classifications into a single csv file using *merge_to_csv.py*
After classification is done. *merge_to_csv.py* will merge all the arff files and comines them 
into a single csv files. As these are easier to work with in the future.

#### Step 5
##### Create a histogram
Uses script: *create_histogram.R*  
A histogram is created and saved. The histogram creates a histogram of ages and the ratio of 
heart disease in each age group. A line will be drawn with the mean age where people have heart disease.

### Running the pipeline
The config.yaml contains all the paths which can be edited. The Snakefile itself 
doesn't need any changes to run. The config file is layout at follows: 

````yaml
# Folders
workdir:
datadir: data/
resultsdir: results/
histogramPath: rendered/histogram.png

# Inputfile
inputfile: heart_data_small.arff

````

### Starting the pipeline
Starting the pipeline is straight forward.
Navigate to workflow directory and use command: 
```commandline
Snakemake -c[amount of cores]
```

## Directory structure
````text
.
├── config
│   └── config.yaml
├── dag.png
├── data
│   ├── heart_data_large.arff
│   ├── heart_data_medium.arff
│   └── heart_data_small.arff
├── readme.md
├── requirements.txt
├── resources
│   ├── adaboost.model
│   └── classify.jar
├── results
│   ├── classified
│   ├── cleaned
│   ├── merged
│   ├── rendered
│   └── split_arff
└── workflow
    ├── Snakefile
    └── scripts
````

## Contact
* S.J. Bouwman
  * s.j.bouwman@st.hanze.nl
