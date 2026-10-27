function RA_Continuum_ValidateOptions(options, D_o)
    % RA_CONTINUUM_VALIDATEOPTIONS Fail early for unsupported settings.
    assert(options.deltaMin > 0, 'deltaMin must be positive.');
    assert(options.deltaMax > options.deltaMin, 'deltaMax must exceed deltaMin.');
    assert(options.deltaMax < 0.5 * D_o, 'deltaMax must leave a positive oil core.');
    assert(options.nScan >= 20, 'Use at least 20 scan points.');
    assert(options.sigma_oa > 0 && options.sigma_oa <= 2, 'Invalid sigma_oa.');
    assert(options.sigma_aw > 0 && options.sigma_aw <= 2, 'Invalid sigma_aw.');
    assert(any(strcmpi(options.buoyancyModel, {'full', 'oil-only'})), ...
           'buoyancyModel must be full or oil-only.');
    assert(any(strcmpi(options.gasClosure, ...
                       {'closed-shell', 'zero-pressure-gradient'})), ...
           'Unknown gasClosure.');
end
