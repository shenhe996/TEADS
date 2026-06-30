## TEADS — TE Anomaly Detection System

# Introduction

The Transposable Element Anomaly Distribution Detection System (TEADS) proposes a novel genome assembly error detection pipeline. It is grounded in biological principles—specifically, the fact that LTR retrotransposons in plants tend to be enriched near chromosome centromeres—to identify severe assembly errors that are otherwise difficult to detect through manual curation or standard sequencing QC. The method targets large-scale structural assembly errors and is particularly suited for assembly error detection in the absence of a reference genome. It also supports batch detection of assembly errors across multiple genomes.

This program provides intuitive visualization of results. Users can also train their own models using custom-annotated training sets by adjusting training parameters, enabling fully self-supervised training and detection.

# Overview

TEADS employs a hybrid detection strategy that integrates deep learning-based reconstruction error analysis with statistical pattern recognition. The system classifies results into three tiers: normal, suspicious, and anomalous. Through built-in visualization tools, it generates TE density distribution plots and detailed statistical reports.

![](https://cdn.nlark.com/yuque/0/2026/png/54929366/1782720796303-45948eaa-2bba-449d-acf3-eba0428184f0.png)

# Installation

## Install the repository and dependencies

```plain
git clone https://github.com/shenhe996/TEADS.git
git clone https://github.com/shenhe996/te_analysis.py
```

It is recommended to create a dedicated conda environment:

```plain
# Create an isolated environment for TEADS
conda create -n TEADS python=3.10
```

Install the following dependencies via conda or pip:

| Package | Required Version | Recommended Version |
| --- | --- | --- |
| numpy | ≥1.21.0 | 2.2.6 |
| pandas | ≥1.3.0 | 2.3.1 |
| matplotlib | ≥3.4.0 | 3.10.3 |
| seaborn | ≥0.11.0 | 0.13.2 |
| scikit-learn | ≥1.0.0 | 1.7.1 |
| torch | ≥2.0.0 | 2.6.0 |

## Install the annotation software

```plain
# Setup
conda create -n EDTA
conda activate EDTA
# Install Mamba inside the EDTA environment
conda install -c conda-forge mamba
mamba install -c conda-forge -c bioconda edta
# Verify installation
EDTA.pl --help
# Verify installation
EDTA.pl --help
```

# Usage

## 1. TE Annotation and Element Density Calculation

Whole-genome TE annotation is obtained using either EDTA or LTR_FINDER combined with LTR_retriever. TE density distributions are then computed using Python 3.12.11. Each chromosome or scaffold is divided into fixed-size windows—typically 100 windows per scaffold, each spanning 3 Mb. For every window, the density of distinct TE categories is calculated, including LTR retrotransposons, LINE elements, SINE elements, DNA transposons, and unclassified TEs.

```plain
# Direct annotation with EDTA
EDTA.pl --genome genome.fa --overwrite 1 --sensitive 1 --anno 1 --threads 20
# Annotate LTR elements only
ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 out_JBAT.FINAL.fa > sunset.finder.scn &
# Faster multi-threaded annotation with LTR_FINDER_parallel
LTR_FINDER_parallel -seq genome_file -threads 20 -harvest_out

# If using ltr_finder output:
LTR_retriever -genome TAIR10.fa -infinder TAIR10.finder.scn -threads 20
# If using LTR_FINDER_parallel output:
LTR_retriever -genome genome_file -inharvest genome_base.finder.combine.scn -threads 20

# After obtaining the .out.gff annotation file, use the te_analysis script to
# compute a summary text file by specifying the appropriate input path.
# For batch processing, files can first be renamed (rename.sh) and then
# processed by ID in bulk using ltranalysis.py, which calculates the
# distribution density of different TE categories along each chromosome.
```

This produces a summary text file.

## Input Data Format

```plain
Chr	From	To	LTR	LINE	SINE	DNA	unknown	TE
scaffold_1	0	100000	0.123	0.045	0.012	0.067	0.003	0.250
scaffold_1	100000	200000	0.134	0.039	0.015	0.071	0.002	0.261
scaffold_1	200000	300000	0.128	0.042	0.013	0.069	0.004	0.256
...
scaffold_2	0	100000	0.089	0.056	0.018	0.045	0.001	0.209
...
```

## 2. Running TEADS

TEADS is run via the Linux command line. Usage and parameters are as follows:

```bash
python TEADS [options]
```

### Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `--mode` | `{train, test, both}` | — | Running mode: train a model, test with an existing model, or both |
| `--train-dir` | `path` | — | Directory containing training data |
| `--test-dir` | `path` | — | Directory containing testing data |
| `--min-length` | `int` | `10000000` | Minimum scaffold length (bp) to be included |
| `--batch-size` | `int` | `16` | Number of samples per training batch |
| `--epochs` | `int` | `100` | Number of training epochs |
| `--lr` | `float` | `0.0001` | Learning rate for model optimization |
| `--z-threshold` | `float` | `4.5` | Z-score threshold for outlier detection |
| `--cnn-doubt-percentile` | `int` | `85` | Percentile threshold for CNN doubt classification |
| `--cnn-anomaly-percentile` | `int` | `90` | Percentile threshold for CNN anomaly classification |

### Examples

```bash
# Train mode
python TEADS --mode train --train-dir /path/to/data/

# Test mode
python TEADS --mode test --test-dir /path/to/data/

# Train and test
python TEADS --mode both --train-dir /path/to/train/ --test-dir /path/to/test/
```

## Output Files

### Generated Files

After running the detection, you'll find:

```plain
te_anomaly_results.csv        # Complete results for all scaffolds
anomalous_scaffolds.csv       # Only anomalous scaffolds
training_curve.png            # Training/validation loss curve (train mode)
best_anomaly_model.pth        # Trained model file
anomaly_all_scaffolds.pdf     # Visualization of LTR distribution anomalies
```

### Result Columns

**te_anomaly_results.csv:**

| Column | Description | Values |
| --- | --- | --- |
| `Scaffold` | Scaffold name | scaffold_1 |
| `CNN_Score` | Reconstruction error from autoencoder | 0.0 - 1.0+ |
| `Pattern_Score` | Distribution pattern anomaly score | 0.0 - 1.0 |
| `Pattern_Type` | Detected pattern type | valley_pattern_anomaly, normal_peak_pattern, etc. |
| `Pattern_Status` | Pattern-based classification | normal, doubt, anomaly |
| `CNN_Status` | CNN-based classification | normal, doubt, anomaly |
| `Has_Outliers` | Whether scaffold has outlier windows | True/False |
| `Final_Status` | Combined final classification | **normal**, **doubt**, **anomaly** |
| `Reason` | Explanation for classification | Contains 3 delta-TE outliers |
| `Source_File` | Original input file | species1_TE_distribution.txt |

### Interpreting the Results

In principle, both "anomaly" and "doubt" classifications warrant closer attention from the user. If a scaffold is flagged with the reason "Both pattern and CNN detected anomaly," it generally corresponds to what we consider a true LTR misdistribution pattern—indicating the presence of a large-scale assembly error on that chromosome, unless the underlying distribution genuinely conforms to such a pattern by nature. Because LTR distributions are inherently diverse, a high CNN score alone or a high pattern score alone does not rule out the possibility of a false positive. However, if manual inspection confirms that the LTR distribution pattern matches a known error type, this strongly indicates that the chromosome does contain an assembly error. Remaining ambiguous cases should be evaluated by integrating alignment against a reference genome and Hi-C signal. The presence of discrete outlier points strongly suggests a breakpoint or an anomalously distributed segment; users may also adjust the detection threshold themselves (default: 4.5).

# Contact

If you have any questions, please feel free to reach out: 2024302010202@webmail.hzau.edu.cn
