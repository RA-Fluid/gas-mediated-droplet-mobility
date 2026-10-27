function index = RA_Continuum_SelectRoot(rootsFound, states, selection)
    % RA_CONTINUUM_SELECTROOT Apply the requested multiple-root rule.
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
