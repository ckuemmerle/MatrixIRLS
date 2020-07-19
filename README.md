# MatrixIRLS
*Matrix Iteratively Reweighted Least Squares* (`MatrixIRLS`) for low-rank matrix completion.

This repository contains a MATLAB implementation of the algorithm  `MatrixIRLS` described in the paper [Escaping Saddle Points in Ill-Conditioned Matrix Completion with a Scalable Second Order Method](https://drive.google.com/file/d/1s-ivhFNLEMe_tSgqUNd-oHD5HCEPiyMF/view) presented at the [ICML 2020 Workshop: _"Beyond First Order Methods in ML Systems"_](https://sites.google.com/view/optml-icml2020/home).

`MatrixIRLS`, as presented in the [paper above]((https://drive.google.com/file/d/1s-ivhFNLEMe_tSgqUNd-oHD5HCEPiyMF/view)),  minimizes the sum of logarithms of the singular values of a matrix subject to a entry-wise data constraint, using Iteratively Reweighted Least Squares (IRLS) steps based on an _optimal weight operator_ combined with a suitable smoothing strategy for the objective.

The implementation uses an adaptation of `bksvd` or [Block Krylov Singular Value Decompostion](https://github.com/cpmusco/bksvd) by Cameron Musco and Christopher Musco to compute singular values and vectors of matrices in a "low-rank + sparse" format. 

The repository contains also a collection of **_reference algorithms_ for low-rank matrix completion**, see list below. In the experimental section of the [paper]((https://drive.google.com/file/d/1s-ivhFNLEMe_tSgqUNd-oHD5HCEPiyMF/view) ), `MatrixIRLS` is compared to these algorithms in terms of _data-efficiency_ (performance for few provided entries) and _scalability_. We provide the implementations of the reference algorithms for the user's convenience in the folder. These implementations all contain only minor modifications of the authors' original code (to allow for timing experiments). Please refer to the respective research papers and original implementations for a description of standard parameter choices.

## Citation
If you refer to the paper or the code, please cite it as:
```
@inproceedings{kuemmerleverdun2020,
  title={Escaping Saddle Points in Ill-Conditioned Matrix Completion with a Scalable Second Order Method},
  author={K{\"u}mmerle, Christian and Verdun, Claudio M.},
  booktitle={Workshop on Beyond First Order Methods in ML Systems at the $37^{th}$ International Conference on Machine Learning},
  year={2020}
}
```
## Installation
* Open MATLAB and run `setup_MatrixIRLS` or, alternatively, add subfolders manually to path. 

The main implementation can be found in the folder `algorithms/MatrixIRLS`, reference algorithms can be found in subfolders of `algorithms`. Some of the algorithms use methods for sparse access of factorized matrices implemented in C compiled as `.mex` files, and other auxiliary files contained in the `tools` folder.
## Examples and Usage
* `demo_MatrixIRLS`: This is a minimal example of how to use the main implementation  `MatrixIRLS.m`.

In the folder `examples/`, you find the following example scripts:
* `example_compare_MC_algos_badlycond`:
Compares several algorihms for low-rank matrix completion for the completion of badly conditioned matrix with condition number $\kappa=:\frac{\sigma_1(\textbf{X}_0)}{\sigma_r(\textbf{X}_0)} = 10^6$.  For this instance, it calculates reconstructions of the algorithms, and visualizes the error decay, both per iteration and with respect to algorithmic runtime.
* `example_compare_MC_algos_wellcond`:
As above, but with a well-conditioned matrix with $\kappa = 2$.
* `example_compare_HMIRLS`: 
Compares `MatrixIRLS` with `HM-IRLS` of the paper [1] for a medium-size problem. The algorithm of  [1] is somewhat similar, but `MatrixIRLS` exhibits improved scalability by orders of magnitude.
* `example_compare_IRLS_various`: 
Compares the iteratively reweighted least squares algorithms `MatrixIRLS`, a similar fast implementation of IRLS
with suboptimal weight operator `'arithmetic'`, the algorithm `HM-IRLS` of [1], and the algorithm `IRLS-col` of [3] based on column reweighting.

In the folder `ICML2020/`, you find the scripts reproducing the experiments of our [ICML Workshop paper](https://sites.google.com/view/optml-icml2020/home).

## List of algorithms
We gratefully acknowledge the authors of the following matrix completion algorithms. For re-use of the algorithms the subfolders of `algorithms/` with the exception of `MatrixIRLS`, please refer to the provided links and contact the authors if the respective terms of usage are unclear.

* `ASD` (_Alternating Steepest Descent_) and `ScaledASD` (_Scaled Alternating Steepest Descent_) by Jared Tanner and Ke Wei, [available](http://www.sdspeople.fudan.edu.cn/weike/code/mc20140528.tar) at [Ke Wei's website](http://www.sdspeople.fudan.edu.cn/weike/publications.html).
* `NIHT` (Normalized Iterative Hard Thresholding) and `CGIHT` (Conjugate Gradient Iterative Hard Thresholding) by Jeffrey Blanchard, Jared Tanner and Ke Wei, [available](http://www.sdspeople.fudan.edu.cn/weike/code/mc20140528.tar) at [Ke Wei's website](http://www.sdspeople.fudan.edu.cn/weike/publications.html).
* `LMaFit` (Low-rank Matrix Fitting algorithm) by Zaiwen Wen, Wotao Yin, and Yin Zhang, available [here](http://lmafit.blogs.rice.edu).
* `LRGeomCG` (Low-rank matrix completion by Geometric CG) by Bart Vandereycken, [available]((http://www.unige.ch/math/vandereycken/matrix_completion.html)) at [Bart Vandereycken's website](http://www.unige.ch/math/vandereycken/research.php).
* `MatrixIRLS`
* `HM-IRLS` (Harmonic Mean Iteratively Reweighted Least Squares) and `AM-IRLS` by Christian Kümmerle and Juliane Sigl, available [here](https://github.com/ckuemmerle/hm_irls).
* `IRLS-col` and `IRLS-row`: IRLS with weight operators that act on the column or row-space, respectively.
* `R2RILS` (Rank 2r Iterative Least Squares) by Jonathan Bauch and Boaz Nadler, available [here](https://github.com/Jonathan-WIS/R2RILS), see also [here](http://www.wisdom.weizmann.ac.il/~nadler/Projects/R2RILS/R2RILS.html). 
* `RTRMC` (Riemannian trust-region method for low-rank matrix completion) by Nicholas Boumal and Pierre-Antoine Absil, available [here](http://web.math.princeton.edu/~nboumal/RTRMC/index.html). The version provided in this repository uses [Manopt 6.0](https://www.manopt.org).
* `ScaledGD` (Scaled Gradient Descent) by Tian Tong, Cong Ma and Yuejie Chi. Available [here](https://github.com/Titan-Tong/ScaledGD).

## Authors
  * Christian Kümmerle
  * Claudio M. Verdun

If you have any problems or questions, please contact us for example via e-mail (kuemmerle@jhu.edu).