function output = HHRTeta(x, eta)
HHR = x.HHR;
taus = x.taus;
[n, k] = size(x.U);
output = eta;
for i = 1 : k
    output = output - taus(i) * HHR(:, i) * (HHR(:, i)' * output);
end
output(1:k,:)=x.s.*output(1:k,:);
end