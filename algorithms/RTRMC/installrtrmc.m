% RTRMC installation script.
% 
% BEFORE your run this script, try to simply execute main and see if
% it works. If it doesn't, then it is probably because the C-Mex codes need
% to be compiled. Launching this script should do the trick. If compilation
% fails, then you probably need to set up your C-compiler correctly for
% Matlab. Type "mex -setup" at the command prompt for help with this.
% From there, instructions for different Matlab versions and operating
% systems are easily reachable.
%
% If you have trouble installing/using this code, feel free to contact the
% authors at: nicolasboumal@gmail.com
%
%
% Please note that you will also need Manopt to use the RTRMC algorithms.
% A version of Manopt should be packaged with this release. If so, simply
% execute the importmanopt script in the corresponding folder. Preferably,
% you can get a more up-to-date version and documentation on
% 
%      www.manopt.org
%
%
% Nicolas Boumal, UCLouvain, October 15, 2014

% Put this flag to 'true' if main failed.
I_launched_main_and_it_failed = true;

if ~I_launched_main_and_it_failed
    error(['Please first try to launch the script main. ' ...
       'If main executes without errors, then '...
       'there is no need to launch installrtrmc. '...
       'If it fails, then edit the flag on line 26 and launch the '...
       'present script again. ']);
end
   
mex -largeArrayDims spmaskmult.c
mex -largeArrayDims setsparseentries.c

% if ~isempty(strfind(computer, 'WIN'))
    % If these don't work, you may try to execute the three other lines...
    mex -largeArrayDims -lmwlapack cholsolvecell.c
    mex -largeArrayDims -lmwlapack -lmwblas spbuildmatrix.c
    mex -largeArrayDims -lmwlapack -lmwblas buildAchol.c
% else
%     % ... right here. The three above have been found to work on at least
%     % one Linux computer, although it seems that the three lines
%     % right here should be the correct way to do it on Linux.
%     mex -largeArrayDims -llapack cholsolvecell.c
%     mex -largeArrayDims -llapack -lblas spbuildmatrix.c
%     mex -largeArrayDims -llapack -lblas buildAchol.c
% end


if exist('manoptsolve', 'file') ~= 2
    cd 'Manopt_2.0';
    importmanopt;
    cd ..;
    fprintf(['Using Manopt 2.0. This might not be ' ...
             'the most up-to-date version of Manopt. Please consider going to ' ...
             'http://www.manopt.org to obtain the latest version.\n\n' ...
             'Manopt 2.0 was added to your Matlab path.\n\n']);
end
