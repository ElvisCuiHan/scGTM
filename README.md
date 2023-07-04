![scGTM Logo](https://github.com/ElvisCuiHan/scKGAM/blob/main/Figures/scGTM.png?width="400")
---


[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg )](https://github.com/ElvisCuiHan/scKGAM/blob/main/LICENSE.md)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

# scGTM: Single-cell generalized trend model

**scGTM** (orignally named as scKGAM) is the abbreviation for *Single-cell Gene Expression Generalized Trend Model*. This is a Python package for modeling the statistical relationship between pseudotime and gene expression data. The paper is available at [bioRXiv](https://www.biorxiv.org/content/10.1101/2021.11.25.470059v1).

It is intended for bioinformatic scientists, applied statisticians, and students who prefer using Metaheuristic algorithms in solving their own bioinformatic optimization problems. scGKM is able to provide various marginal gene distributions with interpretable regression functions. Check out more features!

* **Free software:** MIT license
* **Python versions:** 3.6 and above

## Installation

To install the bleeding-edge version of scGTM, clone this repo:

```shell
$ git clone -b git@github.com:ElvisCuiHan/scGTM.git
```
and then run

```shell
$ cd scGTM
$ python run_scGTM.py --model.iter 100 --model.marginal 'ZIP' --model.save_dir "your/path/to/save" --data.dir "your/path/file.csv" --gene.start 3 --gene.end 4
```

## Usage

scGTM provides a high-level implementation of various marginal distributions including Poisson, negative binomial (NB), zero-inflated Poisson (ZIP) and zero-inflatd negative binomial (ZINB). Further, it utilizes particle swarm optimization algorithm in the package ***pyswarms*** to optimize the objective function. Thus, it aims to be user-friendly and customizable.

The data should be a cell-by-gene matrix where the first column corresponding to the pseudotime:

```math
Index Pseudotime Gene1 Gene2 ...
1.    t1         y11   y12   ...
2.    t2         y21   y22   ... 
3.    t3         y31   y32   ...
4.    t4         y41   y42   ...
```
A typical data structure will be of the following form:

<img src="https://github.com/ElvisCuiHan/scKGAM/blob/main/Figures/data.png" width="700" />

### All-in-one function

Suppose we want to regress Gene 1 on pseudotime using the scGTM, simply we run the `run_scGTM` file in shell:

```shell
python run_scGTM.py --model.iter {# of iterations} --model.marginal 'ZIP' --model.save_dir "your/path/to/save" --data.dir "your/path/file.csv" --gene.start {START INDEX} --gene.end {END INDEX} 
```

and we can replace `run_scGTM.py` with either `run_scGTM_Hill_Only.py` or `run_scGTM_Valley_Only.py` if we are only interested in one of the two trends.

Using the data in our `demo` folder, the command is:

```shell
python run_scGTM_Valley_Only.py --model.iter 120 --model.marginal 'ZIP' --model.save_dir "Demo/Results/" --data.dir "Demo/simu_nb_scGTM_input.csv" --gene.start 1 --gene.end 60
```

- `gene_index`: The index of gene that we want to model.
- `model.marginal`: The marginal distribution of the gene expression, should be one of `["NB", "ZINB", "Poisson", "ZIP"]`.
- `model.iter`: Number of iterations run by PSO, usually 150 suffices.
- `model.save_dir`: The directory to save our results.
- `data.dir`: The path to our data file.
- `gene.start`: Index of the first gene to fit.
- `gene.end`: Index of the last gene to fit.

In the `scGTM.py` file (and the other two), we can modify the arguments to let the model outputs user-defined colors.

- `plot_args`: A dictionary with keys *color* and *cmap*. *color* is a 4x1 vector and *cmap* is a string. For example:
```python
plot_args={
             'color': ['red', 'tomato', 'orange', 'violet'],
             'cmap': 'Blues',
         }
```

If one wants to estimate many genes with different marginals, we can first change the data directory in the function `parallel` and then use the command in terminal:

```bash
python run_scGTM_Hill_Only.py  --gene.start {START INDEX} --gene.end {END INDEX} --model.marginal "NB" --model.save_dir "YourTargetPath" --model.iter 150
```

Note the data should be in *.csv* format. The **main** function will return a *.json* file and *.png* figure. 

### Example

The following figure has shown a typical output by the `main` function in `scGTM.py`.

- *Red line*: fitted log mean expression (log(tau_c) in the paper). 
- *Blue line*: Red line minus -log(1-p_c) so that the zero-inflation part is removed from expectation.
- *Orange vertical line*: Estimated t0, i.e., the turning point of the model.
- *Purple line*: fitted zero-inflation parameter, for details, see paper.
- *Scatters/Points*: observed log expression value (log(y+1)).


<img src="https://github.com/ElvisCuiHan/scKGAM/blob/main/Figures/100ZIP.png" width="700" />

The confidence intervals of `{t0, k1, k2, mu}` are saved in a `.json` file in the same directory.
