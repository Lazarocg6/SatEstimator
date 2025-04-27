function [time,DTtime] = initTimes(inst,duracion,precision,b4andafter)
%INITTIMES Summary of this function goes here
%   Detailed explanation goes here

temp = datenum(inst);

if b4andafter
    time = (temp + (-(duracion/2880):(precision/1440):(duracion/2880)));
else
    time = (temp + (0:(precision/1440):(duracion/1440)));
end

DTtime = datetime(time,'ConvertFrom','datenum','TimeZone','Local');

end

