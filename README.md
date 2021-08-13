![scKGAM Logo](https://github.com/ElvisCuiHan/scKGAM/blob/main/Figures/scKGAM.png?width="400")
---


[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg )](https://github.com/ElvisCuiHan/scKGAM/blob/main/LICENSE.md)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

# scKGAM: Single-cell gene expression kinetics generalized additive model

**scKGAM** is the abbreviation for *Single-cell Gene Expression Kinetics Generalized Additive Model*. This is a Python package for modeling the statistical relationship between pseudotime and gene expression data.

It is intended for bioinformatic scientists, applied statisticians, and students who prefer using Metaheuristic algorithms in solving their own bioinformatic optimization problems. scKGAM is able to provide various marginal gene distributions with interpretable regression functions. Check out more features!

* **Free software:** MIT license
* **Documentation:** https://test.pypi.org/project/scKGAM/1.0/.
* **Python versions:** 3.6 and above

## Installation

To install scKGAM, run this command in your terminal:

```shell
$ pip install -i https://test.pypi.org/simple/ scKGAM==1.0
```

This is the preferred method to install scKGAM. In case you want to install the bleeding-edge version, clone this repo:

```shell
$ git clone -b development https://github.com/ElvisCuiHan/scKGAM.git
```
and then run

```shell
$ cd scKGAM
$ python setup.py install
```

## Usage

scKGAM provides a high-level implementation of various marginal distributions including Poisson, negative binomial (NB), zero-inflated Poisson (ZIP) and zero-inflatd negative binomial (ZINB). Further, it utilizes particle swarm optimization algorithm in the package ***pyswarms*** to optimize the objective function. Thus, it aims to be user-friendly and customizable.

You can import scKGAM as any other Python module,

```python
import scKGAM 
import pyswarms as ps
```

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

Suppose we want to regress Gene 1 on pseudotime using the scKGAM, simply we call the `main` function:

```python
 main(gene_index = 1, marginal="ZIP", iter=50, data_dir=YourDataPath, save_dir=YouTargetPath, plot_args={})
```

- `gene_index`: The index of gene that we want to model.
- `marginal`: The marginal distribution of the gene expression, should be one of `["NB", "ZINB", "Poisson", "ZIP"]`.
- `iter`: Number of iterations run by PSO, usually 150 suffices.
- `data_dir`: The path to our data file.
- `save_dir`: The directory to save our results.
- `plot_args`: A dictionary with keys *color* and *cmap*. *color* is a 4x1 vector and *cmap* is a string. For example:
```python
plot_args={
             'color': ['red', 'tomato', 'orange', 'violet'],
             'cmap': 'Blues',
         }
```

If you want to estimate many genes with different marginals, you can first change the data directory in the function `parallel` and then use the following command in terminal:

```bash
python run_scKGAM.py  --gene.start {START INDEX} --gene.end {END INDEX} --model.marginal "NB" --model.save_dir "YourTargetPath" --model.iter 150
```

Note the data should be in *.csv* format. The **main** function will return a *.json* file and *.png* figure. 

### Example

The following figure has shown a typical output by the `main` function in `scKGAM.py`.

- *Red line*: fitted log mean expression (log(tau_c) in the paper). 
- *Blue line*: Red line minus -log(1-p_c) so that the zero-inflation part is removed from expectation.
- *Orange vertical line*: Estimated t0, i.e., the turning point of the model.
- *Purple line*: fitted zero-inflation parameter, for details, see paper.
- *Scatters/Points*: observed log expression value (log(y+1)).


<img src="https://github.com/ElvisCuiHan/scKGAM/blob/main/Figures/100ZIP.png" width="700" />

The confidence intervals of `{t0, k1, k2, mu}` are saved in a `.json` file in the same directory.
