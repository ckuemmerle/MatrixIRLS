########################################################
#                                                      #
#     Thank you for downloading RTRMC / RCGMC.         #
#                                                      #
########################################################

#
#     To get started: read and run main.m
#

#
#     Feedback and questions welcome: nicolasboumal@gmail.com
#

#
#     This is open source software, distributed under BSD license
#     (See LICENSE.TXT)
#


This zip file was released on July 22, 2015. This is version 3.2.
The only difference compared to 3.1 is moving from Manopt 1.07 to 2.0.


---------------------------------------------------------------------------


This is Matlab code to solve the problem of low-rank matrix completion. The
optimization is executed using optimization on manifolds, via the Manopt toolbox:

    www.manopt.org.

For ease of use, version 2.0 of Manopt is bundled in this zip file, but
it is recommended to download the latest version of Manopt on the website.

The best way to get started is to read and run main.m.
It will generate a random instance of a low-rank matrix completion problem
and solve it. Parameters let the user choose between RCTRMC 1, 2, 2p and
RCGMC and RCGMCp.

An important concept is that of the problem structure. All functions which
manipulate a completion problem will take as input a problem structure.
This structure contains all the necessary information to describe an
instance of the problem: the observation mask, the values of the entries,
the desired rank... It also contains a few preprocessing fields. This is
why, to obtain a problem structure, you must call

  buildproblem.m


The following paper describes the base algorithm,
available in rtrmc.m:

% N. Boumal and P.-A. Absil,
% RTRMC: A Riemannian trust-region method for low-rank matrix completion,
% in the proceedings of the Neural Information Processing Systems conference (NIPS), 2011.


The present code is for a revised and extended version of the algorithm,
as described in the journal paper:

% N. Boumal and P.-A. Absil,
% Low-rank matrix completion via preconditioned optimization on the Grassmann manifold,
% Linear Algebra and its Applications, 2015.
% http://www.sciencedirect.com/science/article/pii/S0024379515001342
% Also available on optimization-online.org:
% http://www.optimization-online.org/DB_HTML/2012/09/3619.html


Both papers are present in the zip file this README came with. See below
for BiBTex entries. Please cite either (or both) of these papers if you use
the accompanying Matlab codes in your research.






@incollection{boumal2011rtrmc,
 title={{RTRMC}: A {R}iemannian trust-region method for low-rank matrix completion},
 author={Boumal, N. and Absil, P.-A.},
 booktitle={Advances in Neural Information Processing Systems 24 ({NIPS})},
 editor={J. Shawe-Taylor and R.S. Zemel and P. Bartlett and F.C.N. Pereira and K.Q. Weinberger},
 pages={406--414},
 year = {2011}
}




@Article{boumal2015rtrmc,
  Title                    = {Low-rank matrix completion via preconditioned optimization on the {G}rassmann manifold},
  Author                   = {Boumal, N. and Absil, P.-A.},
  Journal                  = {Linear Algebra and its Applications},
  Year                     = {2015},
  Pages                    = {200--239},
  Volume                   = {475},

  Doi                      = {10.1016/j.laa.2015.02.027},
  Url                      = {http://www.sciencedirect.com/science/article/pii/S0024379515001342}
}
