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
%         for i=1:instancesize
%             median_err(k,j)=median_err(k,j)+err(k,j,i);
%         end
%         error_avg{k,j}=error_avg{k,j}./instancesize;
    end
end
end

