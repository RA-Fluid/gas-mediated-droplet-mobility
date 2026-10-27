function response = RA_BGK_SolvePlanarResponse(Kn, alpha0, alpha1, opts)
    % Solve three canonical linearised-BGK basis problems
    %
    % Reduced equation:
    %   xi*d(phi)/d(eta) = (u-phi+S)/kappa,
    %   kappa = 2*Kn/sqrt(pi).
    %
    % Basis problems:
    %   U0: lower interface translates with unit speed,
    %   U1: upper interface translates with unit speed,
    %   S : unit pressure-source response.

    arguments
        Kn (1, 1) double {mustBePositive}
        alpha0 (1, 1) double {mustBePositive}
        alpha1 (1, 1) double {mustBePositive}
        opts struct
    end
    assert(alpha0 <= 1 && alpha1 <= 1);
    if abs(alpha0 - 1) > 1e-14 || abs(alpha1 - 1) > 1e-14
        error('Partial accommodation requires a separately validated operator.');
    end

    % Calculate the Gaussian-Hermite polynomials roots
    [xi, w] = RA_BGK_GaussHermite(opts.nVelocity);

    w = w / sqrt(pi); % weight coefficients
    eta = linspace(0, 1, opts.nSpace).'; % non-dimensional normal spacing
    dEta = eta(2) - eta(1); % the spacing size
    kappa = 2 * Kn / sqrt(pi); % kinetic rarefaction parameter

    zero = zeros(opts.nSpace, 1);

    [A, ~] = buildOperator(0, 0, 0);

    M = eye(opts.nSpace) - A;
    basis = [1 0 0; 0 1 0; 0 0 1];
    out = cell(3, 1);
    for b = 1:3
        U0 = basis(b, 1);
        U1 = basis(b, 2);
        S = basis(b, 3);
        [~, bvec] = buildOperator(U0, U1, S);
        u = M \ bvec;
        [uCheck, phi] = transportSweep(u, U0, U1, S);
        err = max(abs(u - uCheck)) / max([1, max(abs(u))]);
        m1 = phi * (w .* xi);
        out{b} = struct('u', u, 'phi', phi, 'q', trapz(eta, u), ...
                        'm10', m1(1), 'm11', m1(end), 'converged', ...
                        isfinite(err) && err <= max(100 * eps, opts.iterationTolerance), ...
                        'error', err);
    end
    tauScale = -sqrt(pi) / Kn;
    response = struct('Kn', Kn, 'kappa', kappa, 'alpha0', alpha0, 'alpha1', alpha1, ...
                      'qU0', out{1}.q, 'qU1', out{2}.q, 'qS', out{3}.q, ...
                      'tau0U0', tauScale * out{1}.m10, 'tau0U1', tauScale * out{2}.m10, ...
                      'tau1U0', tauScale * out{1}.m11, 'tau1U1', tauScale * out{2}.m11, ...
                      'tau0S', tauScale * out{3}.m10, 'tau1S', tauScale * out{3}.m11, ...
                      'converged', all(cellfun(@(z)z.converged, out)), ...
                      'maxFixedPointError', max(cellfun(@(z)z.error, out)), ...
                      'operatorConditionNumber', cond(M), 'basis', {out});

    function [Op, bvec] = buildOperator(U0, U1, S)
        [bvec, ~] = transportSweep(zero, U0, U1, S);
        Op = zeros(opts.nSpace);
        for m = 1:opts.nSpace
            e = zeros(opts.nSpace, 1);
            e(m) = 1;
            [Op(:, m), ~] = transportSweep(e, 0, 0, 0);
        end
    end

    function [uNew, phi] = transportSweep(uSource, U0, U1, S)
        phi = zeros(opts.nSpace, numel(xi));
        target = uSource + S;
        for j = 1:numel(xi)
            if xi(j) > 0
                phi(1, j) = U0;
                E = exp(-dEta / (kappa * xi(j)));
                for iy = 2:opts.nSpace
                    src = 0.5 * (target(iy - 1) + target(iy));
                    phi(iy, j) = E * phi(iy - 1, j) + (1 - E) * src;
                end
            else
                phi(end, j) = U1;
                E = exp(-dEta / (kappa * abs(xi(j))));
                for iy = opts.nSpace - 1:-1:1
                    src = 0.5 * (target(iy + 1) + target(iy));
                    phi(iy, j) = E * phi(iy + 1, j) + (1 - E) * src;
                end
            end
        end
        uNew = phi * w;
    end

end
