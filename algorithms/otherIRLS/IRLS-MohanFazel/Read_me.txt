%%%%%%%%%%%%%%%%%%%%%%
LAST UPDATE: 8/28/2012
%%%%%%%%%%%%%%%%%%%%%%



1. Use Test_IRLS to test the IRLS_q code with randomly generated data or
user-input data.

2. Use Test_sIRLS to test the sIRLS_q code for randomly generated data or
user-input data.

3. Problem Set-Up: Before running Test_IRLS or Test_sIRLS, 
   set-up the problem in "Problem_setup.m"

3a. Type: Choose if you are testing synthetic data or your own data.
    Set Type = 1 for synthetic data and Type = 2 for user-input data.
3b. Synthetic Data: Synthetic data can be generated based on
    the problem instances in "Probinstances.m".
    The problem instances are classified as easy or hard based on
    the # samples available and the degrees of freedom in the matrix.
    i) Choose h = 1 for easy problem instances and h = 2 for hard problem
    instances.
    ii) Finally, choose a specific problem instance: An integer between
        1 and 9.
    
    Some notation used in 'Probinstances.m':
    p: #Measurements or #entries of the matrix that is known. Note p < nxn.
    sr: Sampling ratio = p/(nxn).
    r: Rank of the matrix to be completed.
    fr: Degrees of freedom ratio = rx(2n - r)/p. rx(2n - r) is the
        degrees of freedom in a rank r matrix.
   
3c. User-input data:  
     i) Create a input data matrix M.mat which has the following format:
    "The matrix M.mat has 3 columns. The first two columns denote the row-index and column-index.
     The last column has the values of matrix M at the row-column indices specified in the first two columns." 
     The matrix M represents an incomplete matrix that would be completed by the algorithm.
     ii) Specify an estimate for the rank of the true matrix. 

   
4. Algorithm parameters: 
   The algorithm parameters can be set in "Algorithm_parameters.m",

