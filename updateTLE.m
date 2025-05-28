function updateTLE(lastUpdate,ID)
%UPDATETLE Update TLE from NORAD ID sat array if last update was longer than lastUpdate hours ago 
%   Detailed explanation goes here

    folderName = 'TLEs';
    options = weboptions('Timeout', 30);

    % Check if the folder exists, if not, create it
    if ~exist(folderName, 'dir')
        mkdir(folderName); % Create the folder
    end
    
    fullPath = [fullfile(folderName,int2str(ID)) '.txt'];

    if isfile(fullPath)  % Check if the file exists
    
        fid = fopen(fullPath,'rt');
        fechaUpdate = datetime(fgetl(fid));
        fechaActual = datetime("now");
        diff = hours(fechaActual-fechaUpdate);
        fclose(fid);
    
        if diff > lastUpdate
            fid = fopen(fullPath,'wt');
            fprintf(fid,'%s\n',datetime("now"));
            url = ['https://celestrak.org/NORAD/elements/gp.php?CATNR=',num2str(ID),'&FORMAT=2le'];
            data = webread(url,options);  % Read the content from the URL
            lines = splitlines(data);  % Split the content into lines
            fprintf(fid,'%s\n',lines{1});
            fprintf(fid,'%s\n',lines{2});
            fprintf('TLE %i llamado\n',ID)
            fclose(fid);
        end
        
    else
        fid = fopen(fullPath,'wt');
        fprintf(fid,'%s\n',datetime("now"));
        url = ['https://celestrak.org/NORAD/elements/gp.php?CATNR=',num2str(ID),'&FORMAT=2le'];
        data = webread(url,options);  % Read the content from the URL
        lines = splitlines(data);  % Split the content into lines
        fprintf(fid,'%s\n',lines{1});
        fprintf(fid,'%s\n',lines{2});
        fprintf('TLE %i llamado\n',ID)
        fclose(fid);
    end
end

