% utility: matrix_completion(alg_name,m,n,p,r,epsilon)
%
% Make sure all the called algorithms have the same form of Input and Output
% Input for these algorithms:
%   1) [m,n]: size of the matirx; r: rank of the matrix
%   2) Omega: position of the known entries 
%   3) data: values of the known entries
%   4) start: initilization for the algorithms 
%   5) opts: options for stopping criterior and others 
% Ouputs for the algorithms are formatted into two parts
%   1) the matrix returned by the algorithm
%   2) other information all in a Output structure

function matrix_completion(alg_name,m,n,p,r,epsilon)

if (~exist('alg_name','var')||~exist('m','var')||~exist('n','var')||...
    ~exist('p','var')||~exist('r','var')) 
    error('MC:argChk','Not enough inputs')
end

if ~exist('epsilon','var')
  epsilon = 0; % noise_level = 0
end

% set random stream
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));

% set up the problem set
A = zeros(m*n,1);
A(1:p) = ones(p,1);
[~,ind] = sort(randn(m*n,1));
A = A(ind);
A = reshape(A,m,n);

[I, J] = find(A);
Omega = sub2ind([m n], I, J);

L = randn(m,r);
R = randn(r,n);
M = L * R;
noise = randn(p,1);
noise = epsilon*norm(M(Omega))*noise/norm(noise);
data = M(Omega)+noise;

% default_opts is located in Auxiliary folder
% modify it to use different options.
opts = default_opts;

% convert string to function handle
fhandle = str2func(alg_name);

% make_start_x provide different forms of initialization for different
% algorithms; and it is located in Auxiliary folder
start = make_start_x_IHT_ASD(alg_name,m,n,r,Omega,data);

tic; [Mout, Out] = fhandle(m,n,r,Omega,data, start, opts);t=toc;
relerr = norm(Mout-M,'fro')/norm(M,'fro');
relres = Out.itrelres(end);
reschg = Out.reschg;
iter = Out.iter;

fname = [alg_name '_entry_' datestr(now,'yyyymmdd') '.txt'];
fid = fopen(fname,'a');
fprintf(fid, [alg_name, ' m: %d, n: %d, p: %d, r: %d, noise: %g, t: %g, iter: %d, relerr: %g, relres: %g, reschg: %g\n'], m,n,p,r,epsilon,t,iter,relerr,relres,reschg);
fclose(fid);
