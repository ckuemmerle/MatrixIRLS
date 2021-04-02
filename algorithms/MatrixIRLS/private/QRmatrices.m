function x = QRmatrices(x)
[n, k] = size(x.U);
HHR = zeros(n, k);
taus = zeros(k, 1);
s = zeros(k,1);
A = x.U;

for i = 1 : k
    HHR(i : n, i) = A(i : n, i);
    s(i) = -sign(A(i, i));
    HHR(i, i) = HHR(i, i) + sign(A(i, i)) * norm(A(i : n, i));
    taus(i) = 2 / (HHR(i : n, i)' * HHR(i : n, i));
    A = A - taus(i) * HHR(:, i) * (HHR(:, i)' * A);
end
x.HHR = HHR;
x.taus = taus;
x.s = s;
end