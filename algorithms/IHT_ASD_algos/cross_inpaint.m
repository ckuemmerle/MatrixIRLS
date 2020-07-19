function cross_inpaint(alg_name, img)

X = imread(img);

X=double(X)/255;
r = 50;
[U, S, V] = svds(X,r);

clear X;

M = U*S*V';

[m, n] = size(M);

load('xmask.mat')
[I, J] = find(A);
Omega = sub2ind([m n], I, J);
data=M(Omega);

figure(1)
imshow(A.*M)
drawnow

fprintf('begin reovery? Press Enter\n\n')
pause

% convert string to function handle
fhandle = str2func(alg_name);

opts = default_opts;
start = make_start_x_IHT_ASD(alg_name,m,n,r,Omega,data);

fprintf('recovery starts...\n')

tic;[Mout, Out] = fhandle(m,n,r,Omega,data,start,opts); ttol=toc;
relerr = norm(M-Mout,'fro')/norm(M,'fro');
fprintf([alg_name ' relerr = %g, iter = %d, t = %g, reschg = %g\n'], relerr, Out.iter, ttol, Out.reschg);

fprintf('recovery ends...\n\n')


figure(2)
imshow(Mout)
drawnow

%{
title(alg_name,'fontsize',15);
print([func_name,'_boat_cross'],'-dpdf')
print([func_name,'_boat_cross'],'-depsc')
%}

fprintf('close all the figures? Press Enter\n')
pause
close all
