function updateTLE(lastUpdate,sats)
%UPDATETLE Update TLE from NORAD ID sat array if last update was longer than lastUpdate hours ago 
%   Detailed explanation goes here

    folderName = 'TLEs';


    % Check if the folder exists, if not, create it
    if ~exist(folderName, 'dir')
        mkdir(folderName); % Create the folder
    end

    n_sats = length(sats);

    for n_sat = 1:n_sats
    
        fullPath = [fullfile(folderName,int2str(sats(n_sat).NORAD)) '.txt'];
    
        if isfile(fullPath)  % Check if the file exists
        
            fid = fopen(fullPath,'rt');
            fechaUpdate = datetime(fgetl(fid));
            fechaActual = datetime("now");
            diff = hours(fechaActual-fechaUpdate);
            fclose(fid);
        
            if diff > lastUpdate
                fid = fopen(fullPath,'wt');
                fprintf(fid,'%s\n',datetime("now"));
                url = ['https://celestrak.org/NORAD/elements/gp.php?CATNR=',num2str(sats(n_sat).NORAD),'&FORMAT=2le'];
                data = webread(url);  % Read the content from the URL
                lines = splitlines(data);  % Split the content into lines
                fprintf(fid,'%s\n',lines{1});
                fprintf(fid,'%s\n',lines{2});
                fprintf('TLE %i llamado\n',sats(n_sat).NORAD)
                fclose(fid);
            end
            
        else
            fid = fopen(fullPath,'wt');
            fprintf(fid,'%s\n',datetime("now"));
            url = ['https://celestrak.org/NORAD/elements/gp.php?CATNR=',num2str(sats(n_sat).NORAD),'&FORMAT=2le'];
            data = webread(url);  % Read the content from the URL
            lines = splitlines(data);  % Split the content into lines
            fprintf(fid,'%s\n',lines{1});
            fprintf(fid,'%s\n',lines{2});
            fprintf('TLE %i llamado\n',sats(n_sat).NORAD)
            fclose(fid);
        end
    end
end

