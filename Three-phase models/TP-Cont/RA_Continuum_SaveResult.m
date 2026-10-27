function RA_Continuum_SaveResult(result, options)
    % RA_CONTINUUM_SAVERESULT Save result to a MAT file.
    [pathPart, name, ext] = fileparts(options.outputFile);
    if isempty(ext)
        ext = '.mat';
    end
    if isempty(pathPart)
        pathPart = pwd;
    end
    if ~exist(pathPart, 'dir')
        mkdir(pathPart);
    end
    save(fullfile(pathPart, [name ext]), 'result', '-v7.3');
end
