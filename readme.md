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
##### Clean *ARFF* file from missing values and save to *CSV*
Uses script: *clean_to_csv.R*  
The pipeline uses arff files as input these will be converted to csv files as 
this is an easier file format to work with.  


#### Step 2
##### Classify using *classify_jar_wrapper.py*
Uses script: *parse_data.py* & *classify_jar_wrapper.py*  

First the script *parse_data.py* is called. This script will run *classify_jar_wrapper.py* 
for each instance. The data will be saved in a Pandas Dataframe which will then be saved.

Classification is done using a wrapper for a jar file. 
The jar has some compiling errors therefor the option to the output to file is broken.  
A custom Python wrapper was written that runs the jar and processes the output.

#### Step 3
##### Create a histogram
Uses script: *create_histogram.R*  
A histogram is created and saved. The histogram creates a histogram of ages and the ratio of 
heart disease in each age group. A line will be drawn with the mean age where people have heart disease.

### Running the pipeline
The config.yaml contains all the paths which can be edited. The Snakefile itself 
doesn't need any changes to run. The config file is layout at follows: 

````yaml
workdir: ../data/ 
resultdir: ../results/
cleanedFile: heart_data_cleaned.csv
classifiedCSV: classified.csv
histogramPath: histogram.png
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
├── config.yaml
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
│   ├── classified.csv
│   ├── heart_data_cleaned.csv
│   └── histogram.png
└── workflow
    ├── Snakefile
    ├── classify_jar_wrapper.py
    ├── clean_to_csv.R
    ├── create_histogram.R
    ├── dag.png
    ├── data_visualizer.py
    ├── json_test.py
    └── parse_data.py
````

## Contact
* S.J. Bouwman
  * s.j.bouwman@st.hanze.nl
