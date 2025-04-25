function [time,DTtime] = initTimes(inst,duracion,precision)
%INITTIMES Summary of this function goes here
%   Detailed explanation goes here

temp = datenum(inst);

time = (temp + (0:(precision/1440):(duracion/1440)));

DTtime = datetime(time,'ConvertFrom','datenum');

end

