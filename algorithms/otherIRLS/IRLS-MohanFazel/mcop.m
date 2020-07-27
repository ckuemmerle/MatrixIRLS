function [Z] = mcop(X,m,n,alp,beta);

Z = zeros(m,n); %Pre-allocation
for(i = 1:size(alp,1))
    Z(alp(i),beta(i)) = X(alp(i),beta(i));
end;

return;