Code title: "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion"

Implementation of conjugate gradient (CG) and gradient descent (GD) algorithms 
(c) 2012-2013 Bamdev Mishra <b.mishra@ulg.ac.be>

This package contains a MATLAB implementation of the algorithms presented in the report.

B. Mishra and R. Sepulchre,
"R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
Technical report, arXiv:1306.2672, 2013.

This implementation is due to 
Bamdev Mishra <b.mishra@ulg.ac.be>, 2013.

The implementation is a research prototype still in development and is provided AS IS. 
No warranties or guarantees of any kind are given. Do not distribute this
code or use it other than for your own research without permission of the authors.


Feedback is greatly appreciated.


Installation for R3MC:
-------------------------

- Set current directory as your current folder in Matlab or put it in your Matlab path.
- Run "Install_mex.m". You do not need to do this step for subsequent usage. 
- Run "Run_me_first.m" to add folders to the working path. This needs to be done at the starting of each session.
- To check that everything works, run "Test_R3MC.m" at Matlab command prompt
  (you should see some plots at the end).


Algorithms: 
-----------

- _R3MC: This geometry is to related to the technical report,
 "R3MC: A Riemannian three-factor algorithm for low-rank matrix completion",
  B. Mishra and R. Sepulchre,
  Technical report, arXiv:1306.2672, 2013. 

Files:
------
- _R3MC/ConjugateGradient_URV/R3MC.m                        -------- Main Conjugate Gradient file
- _R3MC/ConjugateGradient_URV/egrad2rgrad.m                 -------- Compute the Riemannian gradient
                             /inner_product_urv.m           -------- Compute inner product between two tangent vectors
                             /proj_tangent_space_urv.m      -------- Project a vector in the Euclidean space onto the horizontal space
                             /proj_horizontal_space_urv.m   -------- Project a vector in the Euclidean space onto the horizontal space
                             /vector_transport_urv.m        -------- Vector transport
- _R3MC/Common_URV/coupled_lyap.m                           -------- Solving the coupled Lyapunov equation arising from horizontal projection                   
                  /tangent_space_lyap.m                     -------- Lyapunov equations arising from tangent space projection                   
                  /guess_linesearch_urv.m                   -------- Linearized stepsize search by solving degree 6 polynomial                
                  /guess_linesearch_urv_accel.m             -------- Accelerated stepsize search by solving an approximate degree 2 polynomial                    



Algorithm details:
------------------
- For Gradient schemes: Four updates for \beta are possible;
                        beta_type = 'off': (steepest) Gradient descent,
                        beta_type = 'P-R': Polak-Riebere+,
                        beta_type = 'F-R': Fletcher-Reeves,
                        beta_type = 'H-S': Hestenes-Stiefel+.

- For linesearch: Two updates are possible;
                  linearized linesearch and 
                  adaptive stepsize (proposed in arXiv:1209.0430).


Disclaimer:
-----------

- All contents are written by Bamdev Mishra (b.mishra@ulg.ac.be) except,
  1) "Auxiliary/*" by Bart Vandereycken from his Matlab code package LRGeom

- Mex files used
  1) updateSparse.c is written by Stephen Becker, 11/10/08
  2) partXY.c is supplied with LMaFit



