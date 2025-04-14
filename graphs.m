function graphs(sats,instante,duracion,precision,RX,filter)
    format long g

    %Filter parameters
    elevMinTX = 10; % degrees
    elevMaxTX = 50;
    azMinTX = 90;
    azMaxTX = 270;
    elevMinRX = 10; % degrees
    elevMaxRX = 50;
    azMinRX = 0;
    azMaxRX = 360;

    freq = 143.050e6;% frequency in hertz
    f = 1/298.26; % WGS72 Parameters
    Re = 6378135; % WGS72 Parameters
    lambda = 3e8/freq;
    
    addpath('SGP4_Vectorized')

    % RADAR eq params

    eirp_tx = 1e6; % EIRP power TX [Watts]
    RCSb = 5; % bistatic RCS [meters^2]

    Grx = 2.15; % Gain RX [dB]
    Lsys = 3; % System loses RX [dB]
    Fs = 3; % Noise figure RX [dB]
    Brx = 0.25e6;% BW RX [Hz]
    Tint = 1; % Integration time [s]

    % GRAVES coords

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

    % Custom functions-----------------------------------------------------------------------------
    % Add UTC time function
    function out_txt = addUTC(event_obj,UTCtime)
        % targetID = find(strcmp({sats.Name}, get(event_obj.Target ...
        %     , 'DisplayName')));
        dat = event_obj.Position;
        index = event_obj.DataIndex;
        out_txt = {['X: ', [num2str(dat(1))] ' min'], ...
                   ['Time: ', char(UTCtime(1, index))], ...
                   ['Y: ', [num2str(dat(2))] '']};

    end

    function out_txt = addUTCsubplot(event_obj, xDataLeft, yDataLeft, xDataRight, yDataRight,UTCtime)

        index = event_obj.DataIndex;
        targetObj = get(event_obj.Target, 'Tag');

        if strcmp(targetObj,'Left')
            out_txt = {['X: ', [num2str(event_obj.Position(1))] ' min'], ...
                       ['Time: ', char(UTCtime(1, index))], ...
                       ['Y: ', [num2str(event_obj.Position(2))] '']};
        elseif strcmp(targetObj,'Right')
            out_txt = {['X: ', [num2str(event_obj.Position(1))] 'MHz'], ...
                       ['Time: ', char(UTCtime(1, index))], ...
                       ['Y: ', [num2str(event_obj.Position(2))] ' min']};
        end
    end

    % Custom function to display time on hover (excluding static markers)
    function output_txt = displayTime(event_obj,t, ...
            latTX, lonTX, elTX, latRX, lonRX, elRX)
    
        % Get the clicked point's tag
        targetObj = get(event_obj.Target, 'Tag');
    
        % Case 1: Trajectory points (show corresponding time)
        if strcmp(targetObj, 'trajectory')
            index = event_obj.DataIndex;
            UData = event_obj.Target.UserData;
            output_txt = {['Time: ', [num2str(UData(2,index)) ' min']], ...
                          ['UTC: ', char(t(index))], ...
                          ['Latitude: ', [num2str(event_obj.Position(1))] 'º'], ...
                          ['Longitude: ', [num2str(event_obj.Position(2))] 'º'], ...
                          ['AASL: ', [num2str(UData(1,index) / 10^3, '%.2f') ' Km']]};
    
        % Case 2: Current satellite position (show time = 0)
        elseif strcmp(targetObj, 'current')
            % targetID = find(strcmp({sats.Name}, erase(get(event_obj.Target ...
            % , 'DisplayName'),' position')));
            UData = event_obj.Target.UserData;
            index = UData(3);
            output_txt = {[erase(get(event_obj.Target, 'DisplayName'),' position')], ...
                          ['Time: ', [num2str(UData(2)) ' min']], ...
                          ['UTC: ', char(t(index))], ...
                          ['Latitude: ', [num2str(event_obj.Position(1))] 'º'], ...
                          ['Longitude: ', [num2str(event_obj.Position(2))] 'º'], ...
                          ['AASL: ', [num2str(UData(1) / 10^3, '%.2f') ' Km']]};
    
        % Case 3: Static points (TX, RX) - Show lat/lon but hide time
        elseif strcmp(targetObj, 'TX')
            output_txt = {['Latitude: ', [num2str(latTX)] 'º'], ...
                          ['Longitude: ', [num2str(lonTX) 'º']], ...
                          ['AASL: ', [num2str(elTX)] ' m']};
    
        elseif strcmp(targetObj, 'RX')
            output_txt = {['Latitude: ', [num2str(latRX)] 'º'], ...
                          ['Longitude: ', [num2str(lonRX)] 'º'], ...
                          ['AASL: ', [num2str(elRX)] ' m']};
    
        % Default: Ignore other objects
        else
            output_txt = {};
        end
    end

    fig3 = figure(3); % If not defined here, the plot resets for each iter
    ax = geoaxes;
    hold(ax, 'on') 
    geolimits([35 50], [-14 14])
    geobasemap(ax, 'darkwater');
    % Plot TX and RX positions
    geoplot(ax, latTX, lonTX, 'xr', 'DisplayName', 'TX', 'Tag', 'TX');
    geoplot(ax, latRX, lonRX, 'xr', 'DisplayName', 'RX', 'Tag', 'RX');

    for num = 1:n_sats
        sat = sats(num);
        % SGP4 propagation
        [recef, vecef, rlla_out, vlla_out, tsince] = propagar(sat, ...
            instante, duracion, precision, filenameEOP, f, Re, NaN);
        
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

        %Gradient descent
        noisyDop = noiseGen(f_doppler);
        optimizedDop = gradDescent(sat, instante, duracion, precision, ...
            filenameEOP, f, Re, latTX, lonTX, elTX, latRX, lonRX, ...
            elRX, lambda, noisyDop);
    
        %AzEl
        [azRX,elevRX,slantRangeRX] = ecef2aer(recef(1,:)*10^3,recef(2,:)*10^3, ...
            recef(3,:)*10^3,RX(1),RX(2),RX(3),referenceEllipsoid('wgs72'));
        [azTX,elevTX,slantRangeTX] = ecef2aer(recef(1,:)*10^3,recef(2,:)*10^3, ...
            recef(3,:)*10^3,latTX,lonTX,elTX,referenceEllipsoid('wgs72'));
    
        if filter
            mask = ~((elevTX>elevMinTX) & (elevRX>elevMinRX) & ...
                (elevTX<elevMaxTX) & (elevRX<elevMaxRX) & ...
                (azTX>azMinTX) & (azRX>azMinRX) & (azTX<azMaxTX) & ...
                (azRX<azMaxRX)); 
        else
            mask = ones(1,length(elevRX));
        end

        % RADAR equation

        l_sys = 10^(Lsys/10);
        g_rx = 10^(Grx/10);
        fs = 10^(Fs/10);

        snr = snr_in(eirp_tx,g_rx,lambda,RCSb,R1,R2,l_sys,fs,Tint);
        
        % Figures----------------------------------------------------------------------------------
        leftX = 0;
        rightX = 1100;
        centerX = 550;
        topY = 550;
        botY = 50; 
        
        % Figure 1---------------------------------------------------------------------------------
        fig1 = figure(1);
        hold on
        plot3(recef(1,:),(recef(2,:)),(recef(3,:)),'DisplayName',sat.Name)
        hold off
        grid on
        pbaspect([1 1 1])
        legend('Location', 'best')
        view(45,45)
        set(fig1, 'Position', [leftX, topY, 560, 420]);
        title('Orbit [Km]')
        
        % Figure 2---------------------------------------------------------------------------------
        % Calculate the distance for each satellite and time point
        aasl = rlla(3,:)/10^3;
        aasl(mask) = NaN;

        fig2 = figure(2);
        grid on
        hold on
        plot(tsince, aasl, 'DisplayName',sat.Name) % Plot vectorized data
        legend('Location', 'best')
        hold off
        title('Altitude above sea level')
        ylabel('Magnitude [Km]')
        xlabel('Time [min]')
        xlim([-duracion duracion])
        set(fig2, 'Position', [rightX, topY, 560, 420]);
        % Set up the data cursor mode with a custom function
        dcm = datacursormode(fig2);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));
        
        % Figure 3---------------------------------------------------------------------------------
        colorOrder = lines(n_sats); % Get default MATLAB color scheme
        satColor = colorOrder(num, :);
        % Plot satellite trajectories
        rlla_t = rlla;
        rlla_t(1,mask) = NaN;
        rlla_t(2,mask) = NaN;

        traject = geoplot(ax, rlla_t(1,:), rlla_t(2,:),'DisplayName' ...
                , sat.Name, 'Tag', 'trajectory','Color',satColor);
        traject.UserData = [rlla(3,:);tsince];

        % Plot current satellite positions
        index = find(tsince == 0);
        curr = geoplot(ax, (rlla(1, index)), (rlla(2, index)), ...
            'ok','DisplayName', [sat.Name ' position'], 'Tag', ...
            'current', 'HandleVisibility', 'off');
        curr.UserData = [rlla(3,index); tsince(index);index];

        % hold(ax, 'off')
        legend(ax, 'show', 'Location', 'southeast');

        % Set up the data cursor mode with a custom function
        dcm = datacursormode(fig3);
        set(dcm, 'UpdateFcn', @(obj, event) displayTime(event, t, ...
            latTX, lonTX, elTX, latRX, lonRX, elRX));

        title(ax,'Current and propagated position')
        set(fig3, 'Position', [leftX, botY, 560, 420]);

        % Figure 4---------------------------------------------------------------------------------
        cte = 1000;
        fig4 = figure(4);
        hold on
        bistaticRange(mask) = NaN;
        plot(tsince, bistaticRange / cte, 'DisplayName',sat.Name);
        hold off
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
        % Set up the data cursor mode with a custom function
        dcm = datacursormode(fig4);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));
        set(fig4, 'Position', [rightX, botY, 560, 420]);
        % 
        % Figure 5---------------------------------------------------------------------------------
        vel = zeros(1, length(tsince));

        % Calculate the distance for each satellite and time point
        for i = 1:length(tsince)
            vel(i) = norm(vecef(:,i));
        end

        fig5 = figure(5);
        grid on
        hold on
        vel(mask) = NaN;
        plot(tsince, vel,'DisplayName',sat.Name) % Plot vectorized data

        legend('Location', 'best')
        hold off
        title('Satellite speed')
        ylabel('Magnitude [Km/s]')
        xlabel('Time [min]')
        xlim([-duracion duracion])
        % Set up the data cursor mode with a custom function
        dcm = datacursormode(fig5);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));    
        set(fig5, 'Position', [centerX, topY, 560, 420]);

        % Figure 6---------------------------------------------------------------------------------
        cte = 1000;
        fig6 = figure(6);
        bistaticVelocity(mask) = NaN;
        hold on
        plot(tsince, bistaticVelocity / cte,'DisplayName',sat.Name)
        hold off
        grid on
        legend('Location', 'best')
        title('Bistatic Velocity')
        ylabel('Magnitude [Km/s]')
        xlabel('Time [min]')
        xlim([-duracion duracion])
        % Set up the data cursor mode with a custom function
        dcm = datacursormode(fig6);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));    
        set(fig6, 'Position', [centerX, botY, 560, 420]);

        % Figure 7---------------------------------------------------------------------------------
        cte = 1000;
        fig7 = figure(7);
        numSeries = size(sats, 2);
        colors = lines(numSeries);
        subplot(1,2,1)
        yyaxis left
        divDop = (f_doppler / cte);
        divDop(mask) = NaN;

        hold on
        h = plot(tsince, divDop,'-','Tag','Left','DisplayName',sat.Name);
        hold off
        ylabel('Doppler frequency [Hz]')
        yyaxis right
        ylim([min(min(f_doppler/cte)) max(max(f_doppler/cte))])
        ytickformat('%.4f') 
        ax2 = gca;
        axCol = 'black';
        ax2.YAxis(1).Color = axCol;
        ax2.YAxis(2).Color = axCol;
        ylim((ylim + freq)/10^6);
        ylabel('RX frequency [MHz]')
        set(h, {'Color'}, num2cell(colors(num,:), 2));
        grid on
        legend('Location', 'best')
        sgtitle('Doppler')
        xlabel('Time [min]')
        xlim([-duracion duracion])
        legend('Location', 'best')

        subplot(1,2,2)
        divDop2 = (freq+f_doppler)/10^6;
        divDop2(mask) = NaN;
        hold on
        plot(divDop2,tsince,'-','Tag','Right','DisplayName',sat.Name);
        hold off
        grid on;
        xlabel('Frequencies [MHz]')
        ylabel('Time [min]')
        ylim([-5 5])
        delta = 0.001;
        % xlabel([(freq/1e6)-delta (freq/1e6)+delta])
        xlabel('Frequency [MHz]')
        ytickformat('%.4f') 
        % Set up the data cursor mode with a custom function
        dcm = datacursormode(fig7);
        set(dcm, 'UpdateFcn', @(obj, event) addUTCsubplot(event, tsince, ...
            f_doppler / cte,(freq+(f_doppler/cte))/10^6,tsince,t));    
        set(fig7, 'Position', [centerX-280, topY, 2*560, 420]);

        % Figure 8---------------------------------------------------------------------------------
        fig8 = figure(8);

        yyaxis left
        hold on
        plot(tsince,(azRX(1,:)),'DisplayName',sat.Name)
        hold off
        ylabel('Azimuth [º]')
        yyaxis right
        hold on
        plot(tsince,(elevRX(1,:)),'DisplayName',sat.Name)
        hold off
        ylabel('Elevation [º]')

        title('Position relative to RX')
        xlabel('Time [min]')
        T = {sats(1:num).Name};
        legend(T(repmat(1:length(T),1,2)),'Location', 'best')
        grid on
        set(fig8, 'Position', [leftX, 3*topY/5, 560, 420]);
        dcm = datacursormode(fig8);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));

        % Figure 9---------------------------------------------------------------------------------
        fig9 = figure(9);

        yyaxis left
        hold on
        plot(tsince',(azTX(1,:)),'DisplayName',sat.Name)
        hold off
        ylabel('Azimuth [º]')
        yyaxis right
        hold on
        plot(tsince',(elevTX(1,:)),'DisplayName',sat.Name)
        hold off
        ylabel('Elevation [º]')

        title('Position relative to TX')
        xlabel('Time [min]')
        T = {sats(1:num).Name};
        legend(T(repmat(1:length(T),1,2)),'Location', 'best')
        grid on
        set(fig9, 'Position', [rightX, 3*topY/5, 560, 420]);
        dcm = datacursormode(fig9);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));

        % Figure 10---------------------------------------------------------------------------------
        fig10 = figure(10);
        
        divSlantTX = (slantRangeTX(1,:))/cte;
        divSlantTX(mask) = NaN;
        divSlantRX = (slantRangeRX(1,:))/cte;
        divSlantRX(mask) = NaN;

        hold on
        plot(tsince,divSlantTX, 'DisplayName',[sat.Name ' TX'])
        hold off
        ylabel('Slant Range [Km]')
        hold on
        plot(tsince,divSlantRX, 'DisplayName',[sat.Name ' RX'])
        hold off
        title('Slant Ranges')
        xlabel('Time [min]')
        % T = {sats(1:num).Name};
        % legend(T(repmat(1:length(T),1,2)),'Location', 'best')
        legend('Location','best')
        grid on
        set(fig10, 'Position', [centerX, topY, 560, 420]);
        dcm = datacursormode(fig10);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));
        % Figure 11---------------------------------------------------------------------------------
        fig11 = figure(11);
        snr(mask) = NaN;
        hold on
        plot(tsince,snr, 'DisplayName',sat.Name)
        hold off
        ylabel('SNR [dB]')
        title('Input signal to noise ratio')
        xlabel('Time [min]')
        legend('Location','best')
        grid on
        set(fig11, 'Position', [centerX, botY, 560, 420]);
        dcm = datacursormode(fig11);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));
        % Figure 12---------------------------------------------------------------------------------
        fig12 = figure(12);
        divDop3 = (freq+(noisyDop/cte))/10^6;
        divDop3(mask) = NaN;
        divDop4 = (freq+(optimizedDop/cte))/10^6;
        divDop4(mask) = NaN;
        hold on
        % plot(divDop2,tsince,'DisplayName','OG Doppler');
        plot(divDop3,tsince,'DisplayName','Noisy Doppler');
        plot(divDop4,tsince,'DisplayName','Optimized Doppler');
        hold off
        grid on;
        xlabel('Frequencies [MHz]')
        ylabel('Time [min]')
        ylim([-5 5])
        xlabel('Frequency [MHz]')
        ytickformat('%.4f') 
        legend('Location','best')
        % Set up the data cursor mode with a custom function
        dcm = datacursormode(fig11);
        set(dcm, 'UpdateFcn', @(obj, event) addUTC(event,t));   
        set(fig12, 'Position', [centerX, 3*topY/5, 560, 420]);

    end
end