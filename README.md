# LUAD analysis using DeepProfiler

This repository contains the source code to run the cmVIP analysis in the LUAD
dataset.

## Steps to reproduce

### 1. Install requirements

This folder is a DeepProfiler project. Follow the [Quick
Guide](https://github.com/cytomining/DeepProfiler/tree/c91b9d821a37d90583d19d209be2e53fe3f08d8d#quick-guide)
to install DeepProfiler.

Note: Experiments reported in the paper used the `c91b9d8` commit.

### 2. Download the data

Be aware this script will override any previous data. To download the data run:

```bash
$ utils/download_all.sh
```

### 3. Prepare the data.

Use DeepProfiler to prepare the dataset:

```bash
$ python3 -m deepprofiler --root=./ --config luad.json --gpu 0 prepare
```

`--gpu` option sets the GPU id to use.

### 4. Extract features.

TODO: Add link to `efficientnet-b0_weights_tf_dim_ordering_tf_kernels_autoaugment.h5`

copy your pretrained model `efficientnet-b0_weights_tf_dim_ordering_tf_kernels_autoaugment.h5`
to the `outputs/efn_pretrained/checkpoint` folder.


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

