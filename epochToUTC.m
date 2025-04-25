function [epochDaytime] = epochToUTC(epoch)
%EPOCHTOUTC Summary of this function goes here
%   Detailed explanation goes here

    fraction = epoch-floor(epoch);
    total_hours = fraction * 24;  % Convert to total hours
    epochH = floor(total_hours);   % Get whole hours
    epochM = floor((total_hours - epochH) * 60);  % Get whole minutes
    epochS = floor(((total_hours - epochH) * 60 - epochM) * 60);  % Get seconds
    epochMS = floor((((total_hours - epochH) * 60 - epochM) * 60 - epochS)*1000);

    year = 2000+floor(epoch/1000);

    dateVal = datetime(year, 1, 1) + days(floor(epoch)-(floor(epoch/1000)*1000) - 1); % Convert DOY to date
    monthNum = month(dateVal);  % Extract the month
    dayNum = day(dateVal, 'dayofmonth');

    epochDaytime = datetime(year,monthNum,dayNum,epochH,epochM,epochS,epochMS,'TimeZone', 'Local');

end

