function median_err = get_median_errors(err)
%get_median_errors

if iscell(err)
   err=cell2mat(err); 
end

[nr_algos,parameter_nrs,~]=size(err);
median_err=zeros(nr_algos,parameter_nrs);
for k=1:nr_algos
    for j=1:parameter_nrs
        median_err(k,j)=median(err(k,j,:));
    end
end
end

