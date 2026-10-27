function RA_BGK_SaveResult(result, options)
    % RA_BGK_SAVERESULT Save one model result to the requested MAT file.
    [pathPart, name, extension] = fileparts(options.outputFile);
    if isempty(extension)
        extension = '.mat';
    end
    if isempty(pathPart)
        pathPart = pwd;
    end
    if ~exist(pathPart, 'dir')
        mkdir(pathPart);
    end
    save(fullfile(pathPart, [name extension]), 'result', '-v7.3');
end
