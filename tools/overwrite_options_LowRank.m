function opts = overwrite_options_LowRank(opts,parameters,j,k)

if isempty(parameters)
    % do nothing
else
    if ~isstruct(parameters)
       error('Parameters must be a struct.') 
    end
    para_names=fieldnames(parameters);
    if strcmp(para_names{1},'r')
        opts.r = parameters.(para_names{1})(j);
    elseif strcmp(para_names{1},'lambda')
        opts.lambda = parameters.(para_names{1})(j);
    elseif strcmp(para_names{1},'m')
        opts.m = parameters.(para_names{1})(j);
    elseif strcmp(para_names{1},'rho')
        opts.rho = parameters.(para_names{1})(j);
    end
    if length(para_names) == 2
        if strcmp(para_names{2},'r')
            opts.r = parameters.(para_names{2})(k);
        elseif strcmp(para_names{2},'lambda')
            opts.lambda = parameters.(para_names{2})(k);
        elseif strcmp(para_names{2},'m')
            opts.m = parameters.(para_names{2})(k);
        elseif strcmp(para_names{2},'rho')
            opts.rho = parameters.(para_names{2})(k);
        end
    end
end
end

