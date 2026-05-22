# Sparse-Tensor-Factorization

````markdown
# Codes for Multi-User Non-Linearly Separable Distributed Computing

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

This repository contains the simulation codes used for the numerical results of the paper:

**Multi-User Non-Linearly Separable Distributed Computing**  
Ali Khalesi, Ahmad Tanha, Derya Malak, and Petros Elia  

The paper is available online on arXiv:

```text
https://arxiv.org/abs/2601.16171v2
````

## Overview

This work studies a multi-user distributed computing problem where several users request the evaluation of non-linearly separable functions. The proposed framework uses sparse tensor factorization, multidimensional tiling, and graph-based assignment methods to reduce the required number of servers under computation and communication constraints.

The codes in this repository reproduce the numerical comparisons reported in the paper.

## Contents

The repository includes Python routines for:

* generating admissible exponent tuples,
* constructing multidimensional tiles,
* checking tuple-tile feasibility,
* building the associated assignment graph,
* solving the assignment problem using max-flow,
* computing the proposed achievable number of servers,
* comparing with the default tensor decomposition baseline,
* comparing with the TDC baseline,
* generating IEEE-style histogram plots.

## Requirements

The codes are written in Python and require:

```bash
numpy
matplotlib
```

The implementation also uses standard Python libraries, including:

```bash
math
itertools
collections
```

## How to Run

To run the simulations and generate the plots, use:

```bash
python main.py
```

The script generates comparison plots for different parameter settings, for example:

```python
K = 6, L = 3, P = 6, Lambda = 3, Delta = 6
K = 8, L = 5, P = 5, Lambda = 2, Delta = 4
```

The generated figures are saved in both PDF and EPS formats.

## Citation

If you use this code or find it useful in your research, please cite our paper:

```bibtex
@article{khalesi2026multiuser,
  title={Multi-User Non-Linearly Separable Distributed Computing},
  author={Khalesi, Ali and Tanha, Ahmad and Malak, Derya and Elia, Petros},
  journal={arXiv preprint arXiv:2601.16171},
  year={2026}
}
```

The paper is available at:

```text
https://arxiv.org/abs/2601.16171v2
```

## Authors

Ali Khalesi
Ahmad Tanha
Derya Malak
Petros Elia

## License

This work is licensed under the **Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License**, also known as **CC BY-NC-ND 4.0**.

You are free to share this material, provided that appropriate credit is given, a link to the license is provided, and any changes are indicated.

You may not use the material for commercial purposes.

You may not distribute modified versions of the material.

For more details, see:

```text
https://creativecommons.org/licenses/by-nc-nd/4.0/
```

## Disclaimer

The codes are provided for research and reproducibility purposes. They are intended to accompany the paper and reproduce the numerical results presented therein.

```
```
