function index = RA_BGK_SelectRoot(rootsFound, states, selection)
    % RA_BGK_SELECTROOT Choose one accepted root when multiple roots exist.
    switch lower(selection)
        case 'smallest'
            [~, index] = min(rootsFound);
        case 'largest'
            [~, index] = max(rootsFound);
        case 'minresidual'
            [~, index] = min(abs([states.residual]));
        otherwise
            error('Unknown rootSelection: %s', selection);
    end
end
