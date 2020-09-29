
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
Compares `MatrixIRLS` with `HM-IRLS` of [[KS18]](http://www.jmlr.org/beta/papers/v19/17-244.html) for a medium-size problem. The algorithm of  [[KS18]](http://www.jmlr.org/beta/papers/v19/17-244.html) is somewhat similar, but `MatrixIRLS` exhibits improved scalability by orders of magnitude.
* `example_compare_IRLS_various`: 
Compares the iteratively reweighted least squares algorithms `MatrixIRLS`, a similar fast implementation of IRLS
with suboptimal weight operator `'arithmetic'`, the algorithm `HM-IRLS` of [[KS18]](http://www.jmlr.org/beta/papers/v19/17-244.html), and the algorithm `IRLS-col` of [[FRW11]](https://epubs.siam.org/doi/abs/10.1137/100811404) based on column reweighting.
* `example_compare_MatrixIRLS_IRucLq_sIRLS`:
Compares `MatrixIRLS` with the IRLS variants `IRLS-p` and `sIRLS-p` of [[MF12]](http://www.jmlr.org/beta/papers/v13/mohan12a.html) and `IRucLq` and `tIRucLq` of [[LXY13]](https://epubs.siam.org/doi/abs/10.1137/110840364) (the latter solves an *unconstrained* formulation, minimizing a rank suggorate + data fit term). Illustrates advantage of `MatrixIRLS` in terms of speed and accuracy.

In the folder `ICML2020/`, you find the scripts reproducing the experiments of our [ICML Workshop paper](https://sites.google.com/view/optml-icml2020/home).

## List of algorithms
We gratefully acknowledge the authors of the following matrix completion algorithms. For re-use of the algorithms the subfolders of `algorithms/` with the exception of `MatrixIRLS`, please refer to the provided links and contact the authors if the respective terms of usage are unclear.

* `ASD` (_Alternating Steepest Descent_) and `ScaledASD` (_Scaled Alternating Steepest Descent_) by Jared Tanner and Ke Wei [[TW16]](https://doi.org/10.1016/j.acha.2015.08.003), [available](http://www.sdspeople.fudan.edu.cn/weike/code/mc20140528.tar) at [Ke Wei's website](http://www.sdspeople.fudan.edu.cn/weike/publications.html).
* `NIHT` (_Normalized Iterative Hard Thresholding_ [[TW12]](https://doi.org/10.1137/120876459)) and `CGIHT` (_Conjugate Gradient Iterative Hard Thresholding_ [[BTW15]](https://doi.org/10.1093/imaiai/iav011)) by Jeffrey Blanchard, Jared Tanner and Ke Wei, [available](http://www.sdspeople.fudan.edu.cn/weike/code/mc20140528.tar) at [Ke Wei's website](http://www.sdspeople.fudan.edu.cn/weike/publications.html).
* `LMaFit` (Low-rank Matrix Fitting algorithm [[WYZ12]](https://doi.org/10.1007/s12532-012-0044-1)) by Zaiwen Wen, Wotao Yin, and Yin Zhang, available [here](http://lmafit.blogs.rice.edu).
* `LRGeomCG` (Low-rank matrix completion by Geometric CG [[V13]](https://doi.org/10.1137/110845768)) by Bart Vandereycken, [available]((http://www.unige.ch/math/vandereycken/matrix_completion.html)) at [Bart Vandereycken's website](http://www.unige.ch/math/vandereycken/research.php).
* `MatrixIRLS` (Matrix Iteratively Reweighted Least Squares). This paper/repository.
* `HM-IRLS` (Harmonic Mean Iteratively Reweighted Least Squares [[KS18]](http://www.jmlr.org/beta/papers/v19/17-244.html)) and `AM-IRLS` by Christian Kümmerle and Juliane Sigl, available [here](https://github.com/ckuemmerle/hm_irls).
* `IRLS-p` and `sIRLS-p`: IRLS with weight operator acting on row space only, solving linear systems by gradient descent [[MF12]](http://www.jmlr.org/beta/papers/v13/mohan12a.html), by Karthik Mohan and Maryam Fazel, [available](https://faculty.washington.edu/mfazel/IRLS_final.zip) at [Maryam Fazel's website](https://faculty.washington.edu/mfazel/).
* `IRucLq` and `tIRucLq` ((truncated) Iterative Reweighted unconstrained Lq for low-rank matrix recovery [[LXY13]](https://epubs.siam.org/doi/abs/10.1137/110840364)) by Zaiwen Wen, Wotao Yin and Yin Zhang, [available](https://xu-yangyang.github.io/codes/IRucLq.zip) at [Yangyang Xu's website](https://xu-yangyang.github.io/papers.html). 
* `IRLS-col` and `IRLS-row`: IRLS with weight operators that act on the column or row space, respectively, and thus very similar to algorithms of [[FRW11]](https://epubs.siam.org/doi/abs/10.1137/100811404) and [[MF12]](http://www.jmlr.org/beta/papers/v13/mohan12a.html). Main purpose: illustrate the influence of the choice of weight operator. 
* `R2RILS` (Rank 2r Iterative Least Squares [[BN20]](https://arxiv.org/abs/2002.01849)) by Jonathan Bauch and Boaz Nadler, available [here](https://github.com/Jonathan-WIS/R2RILS), see also [here](http://www.wisdom.weizmann.ac.il/~nadler/Projects/R2RILS/R2RILS.html). 
* `R3MC` (Riemannian three-factor algorithm for low-rank matrix completion [[MS14]](https://doi.org/10.1109/CDC.2014.7039534)) by Bamdev Mishra and Rodolphe Sepulchre, [available](https://dl.dropboxusercontent.com/s/qzxgax0bg3s8oe2/R3MC_17feb_2017.zip) at [Bamdev Mishra's website](https://bamdevmishra.in/codes/r3mc/). We also included `R3MC-rankupd`, a variant of `R3MC` which optimizes on fixed-rank manifolds with increasing rank (see also [[MS14]](https://doi.org/10.1109/CDC.2014.7039534)).
* `RTRMC` (Riemannian trust-region method for low-rank matrix completion [[BA15]](https://doi.org/10.1016/j.laa.2015.02.027)) by Nicholas Boumal and Pierre-Antoine Absil, available [here](http://web.math.princeton.edu/~nboumal/RTRMC/index.html). The version provided in this repository uses [Manopt 6.0](https://www.manopt.org).
* `ScaledGD` (Scaled Gradient Descent [[TMC20]](https://arxiv.org/abs/2005.08898)) by Tian Tong, Cong Ma and Yuejie Chi. Available [here](https://github.com/Titan-Tong/ScaledGD).

## About this repository
##### Main developer: 
* Christian Kümmerle (<kuemmerle@jhu.edu>)
##### Contributors:
* Claudio M. Verdun (<claudio.verdun@tum.de>)

If you have any problems or questions, please contact us for example via e-mail.

## References
 - [[KS18]](http://www.jmlr.org/beta/papers/v19/17-244.html) Christian Kümmerle and Juliane Sigl, [**Harmonic Mean Iteratively Reweighted Least Squares for Low-Rank Matrix Recovery**](http://www.jmlr.org/beta/papers/v19/17-244.html). _J. Mach. Learn. Res._, 19(47):1–49, 2018.
- [[FRW11]](https://epubs.siam.org/doi/abs/10.1137/100811404) Massimo Fornasier, Holger Rauhut and Rachel Ward, [**Low-rank matrix recovery via iteratively reweighted least squares minimization**](https://epubs.siam.org/doi/abs/10.1137/100811404). _SIAM J. Optim._, 21(4):1614–1640, 2011.
- [[MF12]](http://www.jmlr.org/beta/papers/v13/mohan12a.html) Karthik Mohan and Maryam Fazel, [**Iterative reweighted algorithms for matrix rank minimization**](http://www.jmlr.org/beta/papers/v13/mohan12a.html), _J. Mach. Learn. Res._, 13 (1):3441–3473, 2012.
- [[LXY13]](https://epubs.siam.org/doi/abs/10.1137/110840364) Ming-Jun Lai, Yangyang Xu and Wotao Yin, [**Improved iteratively reweighted least squares for unconstrained smoothed $\ell_q$ minimization**](https://epubs.siam.org/doi/abs/10.1137/110840364), _SIAM J. Numer. Anal._, 51(2):927-957, 2013.
- [[TW16]](https://doi.org/10.1016/j.acha.2015.08.003) Jared Tanner and Ke Wei, [**Low rank matrix completion by alternating steepest descent methods**](https://doi.org/10.1016/j.acha.2015.08.003). _Appl. Comput. Harmon. Anal._, 40(2):417–429, 2016.
- [[TW12]](https://doi.org/10.1137/120876459) Jared Tanner and Ke Wei, [**Normalized Iterative Hard Thresholding for Matrix Completion**](https://doi.org/10.1137/120876459)), _SIAM J. Sci. Comput._, 35(5):S104–S125, 2012.
- [[BTW15]](https://doi.org/10.1093/imaiai/iav011) Jeffrey D. Blanchard, Jared Tanner and Ke Wei, [**CGIHT: conjugate gradient iterative hard thresholding for compressed sensing and matrix completion**](https://doi.org/10.1093/imaiai/iav011), _Inf. Inference_, 4(4):289-327, 2015.
- [[WYZ12]](https://doi.org/10.1007/s12532-012-0044-1) Zaiwen Wen, Wotao Yin and Yin Zhang, [**Solving a low-rank factorization model for matrix completion by a nonlinear successive over-relaxation algorithm**](https://doi.org/10.1007/s12532-012-0044-1), _Math. Prog. Comp._ 4(4):333–361, 2012.
-  [[V13]](https://doi.org/10.1137/110845768) Bart Vandereycken, [**Low-Rank Matrix Completion by Riemannian Optimization**](https://doi.org/10.1137/110845768), _SIAM J. Optim._, 23(2):1214-1236, 2013.
- [[BN20]](https://arxiv.org/abs/2002.01849), Jonathan Bauch and Boaz Nadler, [**Rank 2r iterative least squares: efficient recovery of ill-conditioned low rank matrices from few entries**](https://arxiv.org/abs/2002.01849), _arXiv preprint_, arXiv:2002.01849, 2020.
- [[MS14]](https://doi.org/10.1109/CDC.2014.7039534) Bamdev Mishra and Rodolphe Sepulchre, [**R3MC: A Riemannian three-factor algorithm for low-rank matrix completion**](https://doi.org/10.1109/CDC.2014.7039534), In _53rd IEEE Conference on Decision and Control_, 1137-1142. IEEE, 2014.
- [[BA15]](https://doi.org/10.1016/j.laa.2015.02.027) Nicholas Boumal and Pierre-Antoine Absil, [**Low-rank matrix completion via preconditioned optimization on the Grassmann manifold**](https://doi.org/10.1016/j.laa.2015.02.027), _Linear Algebra Appl._, 15(475):200–239, 2015.
- [[TMC20]](https://arxiv.org/abs/2005.08898) Tian Tong, Cong Ma and Yuejie Chi, [**Accelerating Ill-Conditioned Low-Rank Matrix Estimation via Scaled Gradient Descent**](https://arxiv.org/abs/2005.08898), _arXiv preprint_, arXiv:2005.08898, 2020.
