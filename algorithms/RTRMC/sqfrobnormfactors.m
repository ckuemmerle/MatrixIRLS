function val = sqfrobnormfactors(A, B, X, Y)
% FUNCTION VAL = SQFROBNORMFACTORS(A, B, X, Y)
%
% Inputs:
%   A, X: m-by-r matrices
%   B, Y: r-by-n matrices
%
% Output:
%   val = ||AB-XY||^2_F (squared Frobenius norm of the m-by-n matrix AB-XY)
%
% Complexity:
%   O((m+n)*r^2)
%
% Nicolas Boumal, UCLouvain, April 2013.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: rtrmc

    % Thank you to Bart Vandereycken for pointing out this nice linear
    % algebra trick to compute 'val' with good accuracy even when AB is
    % very close to XY! This is both cheap and numerically robust.
    
    % Notice that AB-XY = [A X] * [B' -Y']'.
    % Compute a thin QR factorization of each term:
    [Q1 R1] = qr([A X], 0);    %#ok
    [Q2 R2] = qr([B' -Y'], 0); %#ok
    % Then, by invariance of the Frobenius norm under the orthonormal
    % transformations Q1 and Q2, we simply get:
    val = norm(R1*R2', 'fro')^2;
    % This is accurate up to 7 digits even close to 0, as compared with
    % quadruple precision computations. Thanks, Bart! :)
    
end
