# R2RILS
This repository contains `Python` and `Matlab` implementations for `R2RILS` as described in J. Bauch and B. Nadler (2020), available in preprint [link](https://arxiv.org/abs/2002.01849), as well as simple demos demonstrating the usage of the `Matlab` and `Python` implementations.
## Usage
#### Python
The entry point to run `R2RILS` is a function with the same name which expects the following parameters:
- X: input matrix to complete. This should be a numpy array of dimensions m x n.
- omega: a mask matrix. 1 on observed entries, 0 otherwise.
- rank: the target rank.
- t_max: maximal number of iterations.
- early_stopping_epsilon: early stopping parameter as described in the paper.
- smart_tol: should increase LSQR's accuracy according to the current quality of the objective. Defaults to True.

This method returns X_hat - `R2RILS` estimate for X0.
#### Matlab
The entry point for running `R2RILS` in Matlab is again a function bearing the same name
<br>
function [X_hat U_hat V_hat] = R2RILS(X, omega, rank, t_max)
<br>
and expects the following parameters:
- X: matrix with observed entries in the set omega.
- omega: array of pairs (i,j) indicating which entries are observed.
- rank: the target rank.
- t_max: maximal number of iterations.

This method returns [X_hat U_hat lambda_hat V_hat, observed_RMSE] where:
- X_hat: rank 2r approximation of X0 (note that if `R2RILS` converges than the limiting point is in fact rank r).
- U_hat: matrix of left singular vectors of the best rank r approximation of X_hat.
- lambda_hat: singular values of the best rank r approximation of X0.
- V_hat: matrix of right singular vectors of the best rank r approximation of X_hat.
