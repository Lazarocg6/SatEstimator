function [out_mat] = evitarSaltos(matr)
%EVITARSALTOS Evita saltos en geoplot introduciendo NaN cuando hay saltos  
%   Detailed explanation goes here
    out_mat = matr;
    for n_sat = 1:size(matr,1)
        for i = 1:size(matr,3)
            if i == 1
                continue;
            elseif abs(matr(n_sat,2,i)-matr(n_sat,2,i-1))>180
                out_mat(n_sat,2,i) = NaN;
            end
        end
    end
end

