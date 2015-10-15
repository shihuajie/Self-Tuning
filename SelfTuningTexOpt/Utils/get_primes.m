function [ P ] = get_primes( N )
%GET_PRIMES return the primes up to N
    P = [];
    for n = 2:N
        if is_prime(n)
            P = [P, n]; 
        end
    end
end

function b = is_prime(n)
    for i = 2:floor(sqrt(n))
        if mod(n, i) == 0
            b = 0;
            return
        end
    end
    b = 1;
end