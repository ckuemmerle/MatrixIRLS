function output = HHReta(x, eta)
HHR = x.HHR;
taus = x.taus;
[n, k] = size(x.U);
output = eta;
for i = k : -1 : 1
    output = output - taus(i) * HHR(:, i) * (HHR(:, i)' * output);
end
end