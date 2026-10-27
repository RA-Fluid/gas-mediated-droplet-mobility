function output = RA_Continuum_MergeStruct(base, user)
    % RA_CONTINUUM_MERGESTRUCT Recursively merge user settings into defaults.
    output = base;
    fields = fieldnames(user);
    for k = 1:numel(fields)
        name = fields{k};
        if isfield(base, name) && isstruct(base.(name)) && isstruct(user.(name))
            output.(name) = RA_Continuum_MergeStruct(base.(name), user.(name));
        else
            output.(name) = user.(name);
        end
    end
end
