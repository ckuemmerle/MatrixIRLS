function weight_vec = get_weight_vec(d1,d2,H,dH,increase_antisymmetricweights)
%get_weight_vec This function returns a weight vector "weight_vec" from the
% information of H and dH. In particular, it covers the case of weight operators 
% with symmetric-skewsymmetric splitting if increase_antisymmetricweights
% == 1.
R=size(H{1},1);
d=min(d1,d2);
D=max(d1,d2);
weight_vec =zeros(R*(d1+d2+R),1);

if increase_antisymmetricweights
    weight_vec(1:R^2)=H{1}(:);
    if d1 == D
        weight_vec((R^2+1):(R*(d1+R)))=...
        reshape(kron(dH{1},ones(1,d1)).',[d1*R,1]);
        weight_vec((R*(d1+R)+1):(R*(d2+d1+R)))  =...
        reshape(kron(dH{2},ones(1,d2)),[d2*R,1]);
    else
        weight_vec((R^2+1):(R*(d1+R)))=...
        reshape(kron(dH{2},ones(1,d1)).',[d1*R,1]);
        weight_vec((R*(d1+R)+1):(R*(d2+d1+R)))  =...
        reshape(kron(dH{1},ones(1,d2)),[d2*R,1]);
    end
else
	weight_vec(1:R^2)=H{1}(:);
    weight_vec((R^2+1):(R*(d1+R)))=...
    reshape(kron(dH{1},ones(1,d1)).',[d1*R,1]);
    weight_vec((R*(d1+R)+1):(R*(d2+d1+R)))  =...
    reshape(kron(dH{1},ones(1,d2)),[d2*R,1]);
end
end

