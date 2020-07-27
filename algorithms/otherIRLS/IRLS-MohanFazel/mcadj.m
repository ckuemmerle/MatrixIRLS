%%%ADJOINT OF THE MATRIX COMPLETION OPERATOR%%%%

function [Z] = macdj(y,n,alp,beta);

Z = zeros(n,n);
for(i = 1:size(alp))
    Z(alp(i),beta(i)) = y(i);
end;

return;