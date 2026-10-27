function tableOut = RA_BGK_BuildResponseTable(fixed, kopts)
    % Precompute dimensionless L-BGK response coefficients

    % in case the Knmin is empty, generate an approximate
    if isempty(kopts.KnMin)
        KnMin = max(2e-3, fixed.lambda_mfp / (0.49 * 2 * fixed.R_1));
    else
        KnMin = kopts.KnMin;
    end

    % in case the Knmax is empty, generate an approximate
    if isempty(kopts.KnMax)
        KnMax = max(KnMin * 1.01, fixed.lambda_mfp / 1e-9);
    else
        KnMax = kopts.KnMax;
    end

    % discretise the range of Kn based on the min, max and number
    % logartihmically
    Kn = logspace(log10(KnMin), log10(KnMax), kopts.nKnTable).';

    % allocate the arrays
    coef = nan(numel(Kn), 9);
    converged = false(size(Kn));
    errors = nan(size(Kn));
    conditionNumber = nan(size(Kn));

    % loop over the Kn values
    for k = 1:numel(Kn)

        % main solver for the planar response
        r = RA_BGK_SolvePlanarResponse(Kn(k), kopts.alpha_oa, kopts.alpha_aw, kopts);

        coef(k, :) = [r.qU0, r.qU1, r.qS, r.tau0U0, r.tau0U1, ...
                      r.tau1U0, r.tau1U1, r.tau0S, r.tau1S];
        converged(k) = r.converged;
        errors(k) = r.maxFixedPointError;
        conditionNumber(k) = r.operatorConditionNumber;
        if kopts.verbose
            fprintf('BGK table %d/%d: Kn %.4g, converged %d, error %.3e\n', ...
                    k, numel(Kn), Kn(k), converged(k), errors(k));
        end
    end

    tableOut = struct('Kn', Kn, 'qU0', coef(:, 1), 'qU1', coef(:, 2), ...
                      'qS', coef(:, 3), 'tau0U0', coef(:, 4), 'tau0U1', coef(:, 5), ...
                      'tau1U0', coef(:, 6), 'tau1U1', coef(:, 7), ...
                      'tau0S', coef(:, 8), 'tau1S', coef(:, 9), ...
                      'converged', converged, 'fixedPointError', errors, ...
                      'operatorConditionNumber', conditionNumber, ...
                      'alpha0', kopts.alpha_oa, 'alpha1', kopts.alpha_aw, ...
                      'nVelocity', kopts.nVelocity, 'nSpace', kopts.nSpace);
end
