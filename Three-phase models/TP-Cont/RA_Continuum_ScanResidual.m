function [scan, brackets] = RA_Continuum_ScanResidual(fixed, options)
    % RA_CONTINUUM_SCANRESIDUAL Scan residual and find sign-changing brackets.
    delta = logspace(log10(options.deltaMin), log10(options.deltaMax), options.nScan).';
    template = compact(RA_Continuum_EmptyState());
    scan = repmat(template, options.nScan, 1);
    for k = 1:options.nScan
        scan(k) = compact(RA_Continuum_EvaluateThickness(delta(k), fixed, options));
    end
    valid = [scan.valid].';
    residual = [scan.residual].';
    idx = find(valid(1:end - 1) & valid(2:end) & isfinite(residual(1:end - 1)) & ...
             isfinite(residual(2:end)) & residual(1:end - 1) .* residual(2:end) <= 0);
    brackets = [delta(idx), delta(idx + 1)];
end

function row = compact(s)
    fields = {'delta_a', 'Kn', 'A_o', 'A_w', 'G_amplitude', 'q_relative', ...
            'b_eq', 'FD_Feng', 'FB_selected', 'residual', 'relativeResidual', ...
            'matrixConditionNumber', 'valid', 'message'};
    row = struct();
    for j = 1:numel(fields)
        row.(fields{j}) = s.(fields{j});
    end
end
