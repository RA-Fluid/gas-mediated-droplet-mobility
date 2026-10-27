function kineticTable = ...
    RA_BGK_LoadOrBuildResponseTable(fixed, kineticOptions)
    % RA_BGK_LOADORBUILDRESPONSETABLE Load cached BGK table or build a new one.
    %
    % kineticTable = RA_BGK_LoadOrBuildResponseTable(fixed,kineticOptions)
    %
    % The function:
    %   1. calculates the resolved Kn table limits;
    %   2. creates a signature describing all table-generating parameters;
    %   3. checks for a compatible cached MAT file;
    %   4. loads the table when the signature matches;
    %   5. otherwise builds and saves a new table.
    %
    % The cache is stored in the same folder as this function. If that folder
    % is not writable, the current working folder is used instead.
    %
    % IMPORTANT
    % Increase signature.cacheVersion whenever the kinetic equations,
    % dimensionalisation, boundary conditions, or numerical algorithm change.

    arguments
        fixed struct
        kineticOptions struct
    end

    %% 1. Resolve the actual Kn range

    [KnMin, KnMax] = resolveKnRange(fixed, kineticOptions);

    %% 2. Generate the table signature

    signature = createSignature( ...
                                kineticOptions, KnMin, KnMax);

    %% 3. Determine cache location

    functionFolder = fileparts(mfilename('fullpath'));

    if isempty(functionFolder)
        functionFolder = pwd;
    end

    cacheName = 'RA_BGK_KineticResponseTable_Cache.mat';
    cacheFile = fullfile(functionFolder, cacheName);

    if ~folderIsWritable(functionFolder)
        warning('RA_BGK:CacheFolderNotWritable', ...
                ['The BGK package folder is not writable. ', ...
                     'The response-table cache will be stored in the ', ...
                     'current working folder.']);

        cacheFile = fullfile(pwd, cacheName);
    end

    %% 4. Attempt to load a compatible table

    if isfile(cacheFile)

        cached = load(cacheFile, ...
                      'kineticTable', 'signature');

        tableExists = isfield(cached, 'kineticTable');
        signatureExists = isfield(cached, 'signature');

        if tableExists && signatureExists

            if isequaln(cached.signature, signature)

                kineticTable = cached.kineticTable;

                validateCachedTable( ...
                                    kineticTable, signature);

                % fprintf([ ...
                %     'Loaded cached BGK response table:\n', ...
                %     '  %s\n', ...
                %     '  Kn range: %.6g to %.6g\n', ...
                %     '  Table points: %d\n'], ...
                %     cacheFile, ...
                %     min(kineticTable.Kn), ...
                %     max(kineticTable.Kn), ...
                %     numel(kineticTable.Kn));

                return

            else

                fprintf([ ...
                         'A BGK cache file was found, but its settings ', ...
                         'do not match the current calculation.\n', ...
                         'A new kinetic table will be generated.\n']);

            end

        else

            fprintf([ ...
                     'The BGK cache file is incomplete. ', ...
                     'A new table will be generated.\n']);

        end
    end

    %% 5. Build a new table

    fprintf([ ...
             'Building a new BGK response table.\n', ...
             '  Kn range: %.6g to %.6g\n', ...
             '  Table points: %d\n', ...
             '  Velocity points: %d\n', ...
             '  Spatial points: %d\n'], ...
            KnMin, KnMax, ...
            kineticOptions.nKnTable, ...
            kineticOptions.nVelocity, ...
            kineticOptions.nSpace);

    % Pass the resolved range explicitly so the table and signature use
    % precisely the same Kn bounds.

    buildOptions = kineticOptions;
    buildOptions.KnMin = KnMin;
    buildOptions.KnMax = KnMax;

    kineticTable = ...
        RA_BGK_BuildResponseTable( ...
                                  fixed, buildOptions);

    validateCachedTable( ...
                        kineticTable, signature);

    %% 6. Save the new table

    try

        save(cacheFile, ...
             'kineticTable', 'signature', '-v7.3');

        fprintf([ ...
                 'Saved BGK response table cache:\n', ...
                 '  %s\n'], cacheFile);

    catch ME

        warning('RA_BGK:CacheSaveFailed', ...
                ['The BGK response table was calculated successfully, ', ...
                     'but the cache could not be saved.\nReason: %s'], ...
                ME.message);

    end

end

function [KnMin, KnMax] = ...
    resolveKnRange(fixed, kineticOptions)
    % RESOLVEKNRANGE Reproduce the table builder's automatic Kn limits.

    if isempty(kineticOptions.KnMin)

        KnMin = max( ...
                    2e-3, ...
                    fixed.lambda_mfp / ...
                    (0.49 * 2 * fixed.R_1));

    else

        KnMin = kineticOptions.KnMin;

    end

    if isempty(kineticOptions.KnMax)

        KnMax = max( ...
                    1.01 * KnMin, ...
                    fixed.lambda_mfp / 1e-9);

    else

        KnMax = kineticOptions.KnMax;

    end

    assert(isfinite(KnMin) && KnMin > 0, ...
           'RA_BGK:InvalidKnMin', ...
           'The resolved minimum Knudsen number is invalid.');

    assert(isfinite(KnMax) && KnMax > KnMin, ...
           'RA_BGK:InvalidKnMax', ...
           'The resolved maximum Knudsen number must exceed KnMin.');

end

function signature = ...
    createSignature(options, KnMin, KnMax)
    % CREATESIGNATURE Record every setting affecting table generation.

    signature = struct();

    % Bump this identifier whenever the governing kinetic implementation
    % changes. This prevents an older table from being reused accidentally.

    signature.cacheVersion = ...
        'RA_BGK_direct_Fredholm_v4_2026_09_24';

    signature.KnMin = KnMin;
    signature.KnMax = KnMax;
    signature.nKnTable = options.nKnTable;

    signature.alpha_oa = options.alpha_oa;
    signature.alpha_aw = options.alpha_aw;

    signature.nVelocity = options.nVelocity;
    signature.nSpace = options.nSpace;

    signature.iterationTolerance = ...
        options.iterationTolerance;

    % These are retained because they remain part of the public options,
    % although the direct Fredholm implementation does not use ordinary
    % relaxed source iteration.

    if isfield(options, 'maxIterations')
        signature.maxIterations = ...
            options.maxIterations;
    else
        signature.maxIterations = NaN;
    end

    if isfield(options, 'relaxation')
        signature.relaxation = ...
            options.relaxation;
    else
        signature.relaxation = NaN;
    end

end

function validateCachedTable(table, signature)
    % VALIDATECACHEDTABLE Verify the loaded or newly generated table.

    requiredFields = { ...
                      'Kn', ...
                      'qU0', ...
                      'qU1', ...
                      'qS', ...
                      'tau0U0', ...
                      'tau0U1', ...
                      'tau1U0', ...
                      'tau1U1', ...
                      'tau0S', ...
                      'tau1S', ...
                      'converged'};

    for k = 1:numel(requiredFields)

        assert(isfield(table, requiredFields{k}), ...
               'RA_BGK:InvalidCachedTable', ...
               'The kinetic table is missing field "%s".', ...
               requiredFields{k});

    end

    assert(numel(table.Kn) == signature.nKnTable, ...
           'RA_BGK:InvalidCachedTable', ...
           'The cached table has an unexpected number of Kn points.');

    assert(all(isfinite(table.Kn)), ...
           'RA_BGK:InvalidCachedTable', ...
           'The cached Kn grid contains nonfinite values.');

    assert(all(diff(table.Kn) > 0), ...
           'RA_BGK:InvalidCachedTable', ...
           'The cached Kn grid is not strictly increasing.');

    assert(all(table.converged), ...
           'RA_BGK:InvalidCachedTable', ...
           'One or more cached BGK table entries did not converge.');

    relativeLowerError = abs( ...
                             min(table.Kn) - signature.KnMin) / ...
        max(signature.KnMin, realmin);

    relativeUpperError = abs( ...
                             max(table.Kn) - signature.KnMax) / ...
        max(signature.KnMax, realmin);

    assert(relativeLowerError < 1e-12, ...
           'RA_BGK:InvalidCachedTable', ...
           'The cached lower Kn limit does not match the current settings.');

    assert(relativeUpperError < 1e-12, ...
           'RA_BGK:InvalidCachedTable', ...
           'The cached upper Kn limit does not match the current settings.');

end

function writable = folderIsWritable(folder)
    % FOLDERISWRITABLE Test whether a temporary file can be created.

    writable = false;

    testFile = fullfile( ...
                        folder, ...
                        sprintf('.RA_BGK_write_test_%d.tmp', ...
                randi(1e9)));

    fileID = fopen(testFile, 'w');

    if fileID < 0
        return
    end

    fclose(fileID);

    if isfile(testFile)
        delete(testFile);
    end

    writable = true;

end
