function updateTLE(lastUpdate,sats,filenameTLEs)
%UPDATETLE Update TLE from NORAD ID sat array if last update was longer than lastUpdate hours ago 
%   Detailed explanation goes here

    n_sats = length(sats);

    if isfile(filenameTLEs)  % Check if the file exists
    
        fid = fopen(filenameTLEs,'rt');
        fechaUpdate = datetime(fgetl(fid));
        fechaActual = datetime("now");
        diff = hours(fechaActual-fechaUpdate);
        fclose(fid);
    
        if diff > lastUpdate
            fid = fopen(filenameTLEs,'wt');
            fprintf(fid,'%s\n',datetime("now"));
            for n_sat = 1:n_sats
                url = ['https://celestrak.org/NORAD/elements/gp.php?CATNR=',num2str(sats(n_sat)),'&FORMAT=2le'];
                data = webread(url);  % Read the content from the URL
                lines = splitlines(data);  % Split the content into lines
                fprintf(fid,'%s\n',lines{1});
                fprintf(fid,'%s\n',lines{2});
            end
            disp('TLEs llamados')
            fclose(fid);
        end
        
    else
        fid = fopen(filenameTLEs,'wt');
        fprintf(fid,'%s\n',datetime("now"));
        for n_sat = 1:n_sats
            url = ['https://celestrak.org/NORAD/elements/gp.php?CATNR=',num2str(sats(n_sat)),'&FORMAT=2le'];
            data = webread(url);  % Read the content from the URL
            lines = splitlines(data);  % Split the content into lines
            fprintf(fid,'%s\n',lines{1});
            fprintf(fid,'%s\n',lines{2});
        end
        disp('TLEs llamados')
        fclose(fid);
    end
end

