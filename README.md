# LUAD analysis using DeepProfiler

This repository contains the source code to run the cmVIP analysis in the LUAD
dataset.

## Profiling

### 1. Install requirements

This folder is a [DeepProfiler](https://github.com/cytomining/DeepProfiler)
project.  Experiments reported in the paper used the
[`c91b9d8`](https://github.com/cytomining/DeepProfiler/tree/c91b9d821a37d90583d19d209be2e53fe3f08d8d#quick-guide)
commit.

To install the dependencies, including the DeepProfiler version we used, run:
```bash
$ pip install -r requirements.txt
```

### 2. Download the data

Be aware this script will override any previous data. To download the data run:

```bash
$ utils/download_all.sh
```

### 3. Prepare the data.

 1. Run `extract_locations.py` script to generate location files.

 2. Use DeepProfiler to prepare the dataset:

```bash
$ python3 -m deepprofiler --root=./ --config luad.json --gpu 0 prepare
```

`--gpu` option sets the GPU id to use.

### 4. Extract features.

Use DeepProfiler to extract features:

```bash
$ python3 -m deepprofiler --gpu 0 --exp efn_pretrained --root ./ --config luad.json profile
```

### 5. Create well profiles.

To create the well-based profiles run:

```bash
$ python3 utils/create_profiles.py
```

It will write a `pd.DataFrame` in parquet with profiles.

## VIP analysis

The analysis is split in three notebooks:
 - [1-Expression-VIP.ipynb](1-Expression-VIP.ipynb): Run the baseline analysis using L1000 profiling.
 - [2-Cell-Morphology-VIP.ipynb](2-Cell-Morphology-VIP.ipynb): Run the Cell Morphology VIP method.
 - [3-Aggregation-plots.ipynb](3-Aggregation-plots.ipynb): Create the plots summarizing results.

## Notes about the dataset

From the paper: 
>  An additional 88 constructs are included in the dataset, representing TP53 alleles that inadvertently had double mutations. A comprehensive description of the process for selecting the constructs that were analyzed is presented in Supplementary Figure 2.

We have filtered out these constructs in the **Filter quality control status** section of the [2-Cell-Morphology-VIP.ipynb](https://github.com/broadinstitute/luad-cell-painting/blob/763436df8fc3d862a6e04e82d0847767f2378c2f/2-Cell-Morphology-VIP.ipynb) notebook.
