function output = RA_BGK_MergeStruct(base, user)
    % Recursively merge user options into defaults.
    output = base;
    fields = fieldnames(user);
    for k = 1:numel(fields)

        name = fields{k};
        if isfield(base, name) && isstruct(base.(name)) && isstruct(user.(name))
            output.(name) = RA_BGK_MergeStruct(base.(name), user.(name));
        else
            output.(name) = user.(name);
        end
    end
