function [ C ] = christoffel(D, q)

q = q(:);

n = size(q);


for k = 1:n
    for j = 1:n
        for i = 1:n
            C(k,j) = element(D, q, i, j, k);
        end
    end
end

C = simplify(C, 100);

end

function [ c ] = element( D, q, i, j, k )
%CHRISTOFFEL Summary of this function goes here
%   Detailed explanation goes here
c = 0.5 * ( diff(D(k,j), q(i)) + diff(D(k, i), q(j)) - diff(D(i, j), q(k)) );
end