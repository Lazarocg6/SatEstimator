function [out_mat] = evitarSaltos(matr)
%EVITARSALTOS Evita saltos en geoplot introduciendo NaN cuando hay saltos  
%   Detailed explanation goes here
    out_mat = matr;
    for i = 1:length(matr)
        if i == 1
            continue;
        elseif abs(matr(2,i)-matr(2,i-1))>180
            out_mat(2,i) = NaN;
        end
    end
end

