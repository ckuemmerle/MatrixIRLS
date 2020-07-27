%% -------------------- SET UP THE PROBLEM ------------------- %%


% ------------------------ KARTHIK MOHAN, EE, UW. --------------- %
% ------------------------ LAST UPDATE: 8/28/2012 --------------- %



%% CHOOSE DATA TYPE - SYNTHETIC DATA OR USER-INPUT DATA

type = 1; % CHOOSE 1 for testing on synthetic data.
          % CHOOSE 2 if you have your own data.
          
          
% PROBLEM SET-UP FOR SYNTHETIC DATA

if  type == 1
h = 1; % CHOOSE 1 for EASY problem instances
       % CHOOSE 2 FOR HARD problem instances
       
ins = 1; % CHOOSE PROBLEM INSTANCE: ENTER an integer between 1 and 9
         % SEE "Probinstances.m" for further details.           
      
end;


%% PROBLEM SET-UP FOR USER-INPUT DATA

if type == 2
                         
    % Load USER data
    
    % NOTE:
    % Create a data matrix M.mat which has the following format: 
    % "The matrix M.mat has 3 columns. The first two columns denote the row-index and column-index.
    %  The last column has the values of matrix M at the row-column indices specified in the first two columns." 
    
    existence = exist('M.mat'); % Check for existence of the DATA matrix.
    if existence == 2
    load('M.mat');
    else
    fprintf('\n\n You have not yet added the input matrix M.mat to the current folder.\n Please do so and run the code again. \n See also M_example.mat for an example file. \n');
    return;
    end;
    
    m = size(M,1); % #rows in the matrix to be comlpeted.
    n = size(M,2); % #columns in the matrix to be completed.
    r = min(10,min(m,n)); %An estimate of the rank of the true matrix
      
end;


    
 
 
    
