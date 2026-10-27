function [scan, brackets] = RA_BGK_ScanResidual(fixed, options, evaluator)
    % RA_BGK_SCANRESIDUAL Evaluate the force residual on a logarithmic thickness grid.
    delta = logspace(log10(options.deltaMin), log10(options.deltaMax), options.nScan).';
    template = RA_BGK_ScanRow(RA_BGK_EmptyState());
    scan = repmat(template, options.nScan, 1);
    for k = 1:options.nScan
        scan(k) = RA_BGK_ScanRow(evaluator(delta(k), fixed, options));
    end
    valid = [scan.valid].';
    residual = [scan.residual].';
    idx = find(valid(1:end - 1) & valid(2:end) & isfinite(residual(1:end - 1)) & ...
               isfinite(residual(2:end)) & residual(1:end - 1) .* residual(2:end) <= 0);
    brackets = [delta(idx), delta(idx + 1)];
end
