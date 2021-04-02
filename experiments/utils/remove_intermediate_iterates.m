function outs = remove_intermediate_iterates(outs)

nr_algos = length(outs);
fields = {'X','Xout','Xhist','Yhist'};
for l=1:nr_algos
    for k=1:length(fields)
        if isfield(outs{l},fields{k}) 
            outs{l} = rmfield(outs{l},fields{k});
        end 
    end
end
end

