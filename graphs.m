function graphs(sats,instante,duracion,precision,RX)
    format long g
    
    freq = 143.050e6;% frequency in hertz
    f = 1/298.26; % WGS72 Parameters
    Re = 6378135; % WGS72 Parameters
    rTierra = 6378.0; % Average radius of Earth in km
    lambda = 3e8/freq;
    
    addpath('SGP4_Vectorized')
    
    latGRAVES = 47.34813145826119;
    lonGRAVES = 5.51487314868131;
    elGRAVES = 180; % Average altitude
    
    % Coordinates to be used
    latTX = latGRAVES;
    lonTX = lonGRAVES;
    elTX = elGRAVES;
    latRX = RX(1);
    lonRX = RX(2);
    elRX = RX(3);
    
    n_sats = length(sats);
    
    %EOP.txt filename
    filenameEOP = 'EOP-All.txt';  
    
    updateEOP(24*2, filenameEOP); % Update EOP every 24 hours
    
    updateTLE(6, sats); % Update TLE every 6 hours
    

    % SGP4 propagation
    [recef, vecef, rlla_out, vlla_out, tsince] = propagar(sats, ...
        instante, duracion, precision, filenameEOP, f, Re);
    
    rlla = evitarSaltos(rlla_out); 
    
    % UTC times
    t = instante+minutes(tsince);
    t.Format = 'dd-MMM HH:mm:ss';
    t.TimeZone = 'Local';

    %Bistatic parameters
    [bistaticRange, bistaticVelocity, R1, R2, llaDIST, ecefDIST] = bistaticParams(latTX, ...
        lonTX, elTX, latRX, lonRX, elRX, recef * 1000, vecef * 1000, f, Re);
    
    %Doppler
    f_doppler = (bistaticVelocity*1000)./-lambda;

    % Add UTC time function--------------------------------------------------------------------

    function out_txt = addUTC(event_obj, xData, yData,UTCtime)

        index = event_obj.DataIndex;
        out_txt = {['X: ', [num2str(xData(index))] ' min'], ...
                   ['Time: ', char(UTCtime(1, index))], ...
                   ['Y: ', [num2str(yData(index))] '']};

    end

    function out_txt = addUTCsubplot(event_obj, xDataLeft, yDataLeft, xDataRight, yDataRight,UTCtime)

        index = event_obj.DataIndex;
        targetObj = get(event_obj.Target, 'Tag');

        if strcmp(targetObj,'Left')
            out_txt = {['X: ', [num2str(xDataLeft(index))] ' min'], ...
                       ['Time: ', char(UTCtime(1, index))], ...
                       ['Y: ', [num2str(yDataLeft(index))] '']};
        elseif strcmp(targetObj,'Right')
            out_txt = {['X: ', [num2str(xDataRight(index))] 'MHz'], ...
                       ['Time: ', char(UTCtime(1, index))], ...
                       ['Y: ', [num2str(yDataRight(index))] ' min']};
        end
    end
    
    % Figures----------------------------------------------------------------------------------
    leftX = 0;
    rightX = 1100;
    centerX = 550;
    topY = 550;
    botY = 50; 
    
    % Figure 1---------------------------------------------------------------------------------
    fig1 = figure(1);
    
    plot3(squeeze(recef(:,1,:))',squeeze(recef(:,2,:))',squeeze(recef(:,3,:))')
    grid on
    pbaspect([1 1 1])
    legend({sats.Name}, 'Location', 'best')
    set(fig1, 'Position', [leftX, topY, 560, 420]);
    title('Orbit [Km]')
    
    % Figure 2---------------------------------------------------------------------------------
    dist = zeros(n_sats, length(tsince));
    
    % Calculate the distance for each satellite and time point
    for j = 1:n_sats
        for i = 1:length(tsince)
            dist(j,i) = norm(recef(j,:,i));
        end
    end
    
    fig2 = figure(2);
    grid on
    hold on
    plot(tsince', dist' - rTierra) % Plot vectorized data
    
    legend({sats.Name}, 'Location', 'best')
    hold off
    title('Altitude above sea level')
    ylabel('Magnitude [Km]')
    xlabel('Time [min]')
    xlim([-duracion duracion])
    set(fig2, 'Position', [rightX, topY, 560, 420]);
    % Set up the data cursor mode with a custom function
    dcm = datacursormode(fig2);
    set(dcm, 'UpdateFcn', @(obj, event) addUTC(event, tsince',dist' - rTierra,t));
    
    % Figure 3---------------------------------------------------------------------------------
    lim = length(rlla) - (length(rlla) - length(recef));
    
    fig3 = figure(3);
    ax = geoaxes; % Create Geographic Axes
    hold(ax, 'on') % Enable hold on for the geographic axes
    geolimits([35 50], [-14 14])
    
    % Plot satellite trajectories
    for i = 1:n_sats
        geoplot(ax, squeeze(rlla(i, 1, 1:lim)), squeeze(rlla(i, 2, 1:lim)), ...
            'DisplayName', sats(i).Name, 'Tag', 'trajectory');
    end
    
    % Plot TX and RX positions
    geoplot(ax, latTX, lonTX, 'xr', 'DisplayName', 'TX', 'Tag', 'TX');
    geoplot(ax, latRX, lonRX, 'xr', 'DisplayName', 'RX', 'Tag', 'RX');
    
    % Plot current satellite positions
    for i = 1:n_sats
        index = find(tsince(i, :) == 0);
        geoplot(ax, squeeze(rlla(i, 1, index)), squeeze(rlla(i, 2, index)), ...
            'ok','DisplayName', [sats(i).Name ' position'], 'Tag', ...
            'current', 'HandleVisibility', 'off');
    end
    
    geobasemap(ax, 'darkwater');
    hold off
    legend('Location', 'southeast')
    
    % Set up the data cursor mode with a custom function
    dcm = datacursormode(fig3);
    set(dcm, 'UpdateFcn', @(obj, event) displayTime(event, rlla, tsince, ...
        t, sats, latTX, lonTX, elTX, latRX, lonRX, elRX));
    
    % Custom function to display time on hover (excluding static markers)
    function output_txt = displayTime(event_obj, rlla, tsince,t, sats, ...
            latTX, lonTX, elTX, latRX, lonRX, elRX)
    
        % Get the clicked point's tag
        targetObj = get(event_obj.Target, 'Tag');
        targetID = find(strcmp({sats.Name}, erase(get(event_obj.Target ...
            , 'DisplayName'), ' position')));
    
        % Case 1: Trajectory points (show corresponding time)
        if strcmp(targetObj, 'trajectory')
            index = event_obj.DataIndex;
            output_txt = {['Time: ', [num2str(tsince(targetID, index)) ' min']], ...
                          ['UTC: ', char(t(targetID, index))], ...
                          ['Latitude: ', [num2str(rlla(targetID, 1, index))] 'º'], ...
                          ['Longitude: ', [num2str(rlla(targetID, 2, index))] 'º'], ...
                          ['Elevation: ', [num2str(rlla(targetID, 3, index) / 10^3, '%.2f') ' Km']]};
    
        % Case 2: Current satellite position (show time = 0)
        elseif strcmp(targetObj, 'current')
            index = find(tsince(targetID, :) == 0);
            output_txt = {['Time: ', [num2str(tsince(targetID, index)) ' min']], ...
                          ['UTC: ', char(t(targetID, index))], ...
                          ['Latitude: ', [num2str(rlla(targetID, 1, index))] 'º'], ...
                          ['Longitude: ', [num2str(rlla(targetID, 2, index))] 'º'], ...
                          ['Elevation: ', [num2str(rlla(targetID, 3, index) / 10^3, '%.2f') ' Km']]};
    
        % Case 3: Static points (TX, RX) - Show lat/lon but hide time
        elseif strcmp(targetObj, 'TX')
            output_txt = {['Latitude: ', [num2str(latTX)] 'º'], ...
                          ['Longitude: ', [num2str(lonTX) 'º']], ...
                          ['Elevation: ', [num2str(elTX)] ' m']};
    
        elseif strcmp(targetObj, 'RX')
            output_txt = {['Latitude: ', [num2str(latRX)] 'º'], ...
                          ['Longitude: ', [num2str(lonRX)] 'º'], ...
                          ['Elevation: ', [num2str(elRX)] ' m']};
    
        % Default: Ignore other objects
        else
            output_txt = {};
        end
    end
    
    title('Current and propagated position')
    set(fig3, 'Position', [leftX, botY, 560, 420]);
    
    % Figure 4---------------------------------------------------------------------------------
    cte = 1000;
    fig4 = figure(4);
    plot(tsince', bistaticRange' / cte)
    grid on
    % hold on %Uncomment to show ranges
    % plot(tsince(1,:), R1(1,:) / cte)
    % plot(tsince(1,:), R2(1,:) / cte)
    % hold off
    legend('Location', 'best')
    title('Bistatic range')
    ylabel('Magnitude [Km]')
    xlabel('Time [min]')
    xlim([-duracion duracion])
    legend({sats.Name}, 'Location', 'best')
    % Set up the data cursor mode with a custom function
    dcm = datacursormode(fig4);
    set(dcm, 'UpdateFcn', @(obj, event) addUTC(event, tsince', ...
        bistaticRange'/ cte,t));
    set(fig4, 'Position', [rightX, botY, 560, 420]);
    
    % Figure 5---------------------------------------------------------------------------------
    vel = zeros(n_sats, length(tsince));
    
    % Calculate the distance for each satellite and time point
    for j = 1:n_sats
        for i = 1:length(tsince)
            vel(j,i) = norm(vecef(j,:,i));
        end
    end
    
    fig5 = figure(5);
    grid on
    hold on
    plot(tsince', vel') % Plot vectorized data
    
    legend({sats.Name}, 'Location', 'best')
    hold off
    title('Satellite speed')
    ylabel('Magnitude [Km/s]')
    xlabel('Time [min]')
    xlim([-duracion duracion])
    % Set up the data cursor mode with a custom function
    dcm = datacursormode(fig5);
    set(dcm, 'UpdateFcn', @(obj, event) addUTC(event, tsince',vel',t));    
    set(fig5, 'Position', [centerX, topY, 560, 420]);
    
    % Figure 6---------------------------------------------------------------------------------
    cte = 1000;
    fig6 = figure(6);
    plot(tsince', bistaticVelocity' / cte)
    grid on
    legend('Location', 'best')
    title('Bistatic Velocity')
    ylabel('Magnitude [Km/s]')
    xlabel('Time [min]')
    xlim([-duracion duracion])
    legend({sats.Name}, 'Location', 'best')
    % Set up the data cursor mode with a custom function
    dcm = datacursormode(fig6);
    set(dcm, 'UpdateFcn', @(obj, event) addUTC(event, tsince', ...
        bistaticVelocity'/cte,t));    
    set(fig6, 'Position', [centerX, botY, 560, 420]);
    
    % Figure 7---------------------------------------------------------------------------------
    cte = 1000;
    fig7 = figure(7);
    numSeries = size(f_doppler, 1);
    colors = lines(numSeries);
    subplot(1,2,1)
    yyaxis left
    h = plot(tsince', f_doppler' / cte,'-','Tag','Left');
    ylabel('Doppler frequency [Hz]')
    yyaxis right
    ylim([min(min(f_doppler/cte)) max(max(f_doppler/cte))])
    ytickformat('%.4f') 
    ax = gca;
    axCol = 'black';
    ax.YAxis(1).Color = axCol;
    ax.YAxis(2).Color = axCol;
    ylim((ylim + freq)/10^6);
    ylabel('RX frequency [MHz]')
    set(h, {'Color'}, num2cell(colors, 2));
    grid on
    legend('Location', 'best')
    sgtitle('Doppler')
    xlabel('Time [min]')
    xlim([-duracion duracion])
    legend({sats.Name}, 'Location', 'best')
    
    subplot(1,2,2)
    plot((freq+(f_doppler'/cte))/10^6,tsince','-','Tag','Right');
    grid on;
    xlabel('Frequencies [MHz]')
    ylabel('Time [min]')
    ylim([-10 10])
    ytickformat('%.4f') 
    % Set up the data cursor mode with a custom function
    dcm = datacursormode(fig7);
    set(dcm, 'UpdateFcn', @(obj, event) addUTCsubplot(event, tsince', ...
        f_doppler' / cte,(freq+(f_doppler'/cte))/10^6,tsince',t));    
    set(fig7, 'Position', [centerX-280, 3*topY/5, 2*560, 420]);
end