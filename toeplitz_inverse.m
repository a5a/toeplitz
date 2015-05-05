function [c_inv, logdet] = toeplitz_inverse(m)
%ASA Returns inverse of a Toeplitz matrix and log det c
%   Input is the first column of the Toeplitz matrix
    c = m(:, 1)';
    N = length(c);
    c_inv = zeros(N);
    [v, logdet] = alg1_2(c);
    
    c_inv(1, 1:N) = v(N:-1:1);
    c_inv(1:N, 1) = v(N:-1:1);
    c_inv(N, 1:N) = v(1:N);
    c_inv(1:N, N) = v(1:N);
    
    for i = 2 : floor((N-1)/2) + 1
        for j = i : N - i + 1
            c_inv(i, j) = c_inv(i-1, j-1) + (v(N+1-j) * v(N+1-i) - v(i-1) * v(j-1))/v(N);
%             c_inv(i, j)
            c_inv(j, i) = c_inv(i, j);
            c_inv(N-i+1, N-j+1) = c_inv(i, j);
            c_inv(N-j+1, N-i+1) = c_inv(i, j);
        end
    end
end



function [v, l] = alg1_2(c)
    N = length(c);
    wiggle = c(2:end)'/c(1);
    [z, l] = alg1_3(N-1, wiggle);
    l = l + N * log(c(1));
    v(N) = 1/((1 + wiggle' * z) * c(1));
    v(1:N-1) = v(N) * z(N-1:-1:1);
end


function [z, l] = alg1_3(m, wiggle)
    z = zeros(m, 1);
    z(1) = -wiggle(1);
    beta = 1;
    alpha = -wiggle(1);
    l = 0;
    for i = 1:m-1
        beta = (1-alpha^2) * beta;
        l = l + log(beta);
        alpha = -(wiggle(i+1) + wiggle(i:-1:1)' * z(1:i)) / (beta);
        z(1:i) = z(1:i) + alpha * z(i:-1:1);
        z(i+1) = alpha;
    end
    beta = (1-alpha^2) * beta;
    l = l + log(beta);
end