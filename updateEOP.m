function updateEOP(lastUpdate,filenameEOP)
%UPDATEEOP Update EOP if last update was longer than lastUpdate hours ago 
%   Detailed explanation goes here

    if isfile(filenameEOP)  % Check if the file exists
        fid = fopen(filenameEOP, 'r');  % Open for reading only
        if fid == -1
            error('Error opening the file.');
        else
            version = fgetl(fid);
            fecha = fgetl(fid);
        end
    
        fechaSplitted = split(fecha);
        horaSplitted = split(fechaSplitted{5},':');
        
        months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
              
        % Find the index of the input month in the list
        monthNum = find(strcmpi(fechaSplitted{3}, months));
    
        fechaUpdate = datetime(str2double(fechaSplitted{2}),monthNum ...
            ,str2double(fechaSplitted{4}),str2double(horaSplitted{1}) ...
            ,str2double(horaSplitted{2}),str2double(horaSplitted{3}));
        fechaActual = datetime("now");
    
        if hours(fechaActual-fechaUpdate) > lastUpdate
            websave(filenameEOP,'https://celestrak.org/SpaceData/EOP-All.txt');
            disp('EOP llamado \n')
        end
    
        fclose(fid);
    else
        websave(filenameEOP,'https://celestrak.org/SpaceData/EOP-All.txt');
        disp('EOP llamado \n')
    end 

end

