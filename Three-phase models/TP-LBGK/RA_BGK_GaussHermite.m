function [x, w] = RA_BGK_GaussHermite(n)
    % Nodes and weights for exp(-x^2) Gauss-Hermite quadrature.

    i = (1:n - 1).';

    beta = sqrt(i / 2);

    J = diag(beta, 1) + diag(beta, -1);

    [V, D] = eig(J);

    [x, index] = sort(diag(D));
    V = V(:, index);

    w = sqrt(pi) * (V(1, :).^2).';
