function [ L, a, b ] = labmap( N, mcol )
%LABMAP Create L, a and b color maps

    if nargin < 1
        N = 64;
    end

    Nm1 = N - 1;

    % Center color
    if nargin < 2
        mcol = [0, 0, 0];
    end
    C = repmat(mcol, N, 1);

    % L: black to white
    L = gray(N);
    % midcolor for a and b
    m = 1 + Nm1 / 2;
    m0 = floor(m);
    m1 = m0 + 1;
    if m0 == m1
        m1 = m0 + 1;
    end
    M0 = interp1([1, m], [1, 0], 1:m0);
    M0 = repmat(M0, 3, 1)';
    M1 = interp1([m, N], [0, 1], m1:N);
    M1 = repmat(M1, 3, 1)';
    % a: green to magenta
    green = [0, 1, 0];
    magenta = [1, 0, 1];
    G = repmat(green, N, 1);
    M = repmat(magenta, N, 1);
    
    a = [ M0 .* G(1:m0, :) + (1-M0) .* C(1:m0, :); ...
          M1 .* M(m1:end, :) + (1-M1) .* C(m1:end, :) ];
    % b: blue to yellow
    blue = [0, 0, 1];
    yellow = [1, 1, 0];
    B = repmat(blue, N, 1);
    Y = repmat(yellow, N, 1);
    b = [ M0 .* B(1:m0, :) + (1-M0) .* C(1:m0, :); ...
          M1 .* Y(m1:end, :) + (1-M1) .* C(m1:end, :) ];
end

