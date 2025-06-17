load ../real_data/MatLab_code_Modified_NIF/example/Cascadia_xubo_setup4_f1b1W1_.mat;
load ../real_data/MatLab_code_Modified_NIF/example/GPS_info.mat;
load ../real_data/MatLab_code_Modified_NIF/example/coastline.dat;
%% Plot the slip and slip rates
beginslip1 = 8; 
endslip1 = 80; 
skip1 = 8;
dur1 = 5;
beginslip2 = 6; 
endslip2 = 80; 
skip2 = 10;
dur2 = 10;
xlims = [-125.6 -121.8];
ylims = [43.8 48.9];

slip_diff_sigma_no_up_large_mesh = readmatrix("saved_data/slip_FMOU_large_mesh.csv");
slip_rate_diff_sigma_no_up_large_mesh = readmatrix("saved_data/slip_rate_FMOU_large_mesh.csv");
slip_NIF_large_mesh = readmatrix("saved_data/slip_NIF_large_mesh.csv");
slip_rate_NIF_large_mesh = readmatrix("saved_data/slip_rate_NIF_large_mesh.csv");
slip_NIF_large_mesh_select_d = readmatrix("saved_data/slip_NIF_large_mesh_select_d.csv");
slip_rate_NIF_large_mesh_select_d = readmatrix("saved_data/slip_rate_NIF_large_mesh_select_d.csv");
tremor = readmatrix("../real_data/Cascadia_data/tremor_2011.csv");






dur = 7;
xlims = [-125.6 -121.8];
ylims = [43.8 48.9];


plot_list = [7, 38, 69];
figure
c = 0;
%numplots=ceil((endslip3-beginslip3)/skip3);
numplots=length(plot_list);
h = zeros(1, 3*numplots);

%for i = beginslip3:skip3:endslip3
for i = plot_list
    day1 = stationstruct(i).DOY;
    day2 = stationstruct(i+dur-1).DOY;
    year1=stationstruct(i).year;
    year2=stationstruct(i+dur).year;
    date1 = datetime(2011, 1, 1) + days(day1 - 1);  
    date2 = datetime(2011, 1, 1) + days(day2 - 1); 
    datestr1 = datestr(date1, 'mmm dd');
    datestr2 = datestr(date2, 'mmm dd');
 
    epoch = i;
    %slip rate
     slip_rate = mean(slip_rate_diff_sigma_no_up_large_mesh(:,epoch:(epoch+dur-1)), 2)*100/365;
%     if c == 1
%        slip_rate = mean(slip_rate_NIF_large_mesh(:,epoch:(epoch+dur-1)), 2)*100/365;
%     end
%     if c == 2
%        slip_rate = mean(rate_b(:,epoch:(epoch+dur-1)), 2)*100/365;
%     end

    c=c+1;
   
    h(c) = subplot(1, 3*ceil(numplots), c);
    hold on

    %trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip_rate, 'EdgeAlpha', 0.225);

    validVertices = nd_ll(:,3) <= -20 & nd_ll(:,3) >= -60;
    validElements = all(validVertices(el), 2);
    filteredel = el(validElements, :);
    filteredsliprate = slip_rate(validElements, :);
    trisurf(filteredel,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),filteredsliprate, 'EdgeAlpha', 0.15);

    %caxis([-.25, .25])
    caxis([0, .25]);
    view(2);
    daspect([1/cosd(abs(origin(1))),1,110]);

    hold on;
    tremor_here = tremor(tremor(:,8) >= day1 & tremor(:,8) <= day2 & tremor(:,1) < 48.5, :);
    %plot(tremor_here(:,2), tremor_here(:,1), 'k.', 'MarkerSize', 4.4);
    
    %plot3(tremor_here(:,2), tremor_here(:,1), tremor_here(:,3), '.', 'MarkerSize', 0.01, 'Color', [1, 0.7, 0.85]);
    x = tremor_here(:, 2);
    y = tremor_here(:, 1);
    z = tremor_here(:, 3);
    scatter3(x, y, z, 2, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8, 0, 0.8], 'MarkerFaceAlpha', 0.2);

    %xlabel('Longitude', 'fontsize', 14); ylabel('Latitude', 'fontsize', 14);
   
    %title(['slip-rate, cm/day, days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 12 );
    %title(['days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 10 );
    title([datestr1, ' - ',datestr2, ' ' ], 'fontsize', 10 );

    plot(coastline(:, 1), coastline(:, 2), 'k', 'Linewidth', 0.65)

    ylim(ylims)
    xlim(xlims)
    set(gca, 'fontsize', 10)
    if c > 1
        %set(h(c), 'YTickLabel', [], 'YTick', []);   
        set(h(c), 'YTickLabel', []);   
    else
        xlabel("Longtitude")
        ylabel("Latitude")
    end

    % Reduce the spacing by adjusting the subplot position
    pos = get(h(c), 'Position');
    pos(1) = pos(1) - (c-1) * 0.01; % Shift to the left; adjust the factor as needed
    %pos(2) = pos(2) - 0.01;
    pos(3) = pos(3) + 0.006;         % Increase width; adjust the factor as needed
    set(h(c), 'Position', pos);
end

c = numplots;
for i = plot_list
    day1 = stationstruct(i).DOY;
    day2 = stationstruct(i+dur-1).DOY;
    year1=stationstruct(i).year;
    year2=stationstruct(i+dur).year;
    date1 = datetime(2011, 1, 1) + days(day1 - 1);  
    date2 = datetime(2011, 1, 1) + days(day2 - 1); 
    datestr1 = datestr(date1, 'mmm dd');
    datestr2 = datestr(date2, 'mmm dd');
 
    epoch = i;
    %slip rate
    
    slip_rate = mean(slip_rate_NIF_large_mesh(:,epoch:(epoch+dur-1)), 2)*100/365;

    c=c+1;
   
    h(c) = subplot(1, 3*ceil(numplots), c);
    hold on

    %trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip_rate, 'EdgeAlpha', 0.225);

    validVertices = nd_ll(:,3) <= -20 & nd_ll(:,3) >= -60;
    validElements = all(validVertices(el), 2);
    filteredel = el(validElements, :);
    filteredsliprate = slip_rate(validElements, :);
    trisurf(filteredel,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),filteredsliprate, 'EdgeAlpha', 0.15);

    %caxis([-.25, .25])
    caxis([0, .25]);
    view(2);
    daspect([1/cosd(abs(origin(1))),1,110]);

    hold on;
    tremor_here = tremor(tremor(:,8) >= day1 & tremor(:,8) <= day2 & tremor(:,1) < 48.5, :);
    %plot(tremor_here(:,2), tremor_here(:,1), 'k.', 'MarkerSize', 4.4);
    
    %plot3(tremor_here(:,2), tremor_here(:,1), tremor_here(:,3), '.', 'MarkerSize', 0.01, 'Color', [1, 0.7, 0.85]);
    x = tremor_here(:, 2);
    y = tremor_here(:, 1);
    z = tremor_here(:, 3);
    scatter3(x, y, z, 2, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8, 0.0, 0.8], 'MarkerFaceAlpha', 0.2);

    %xlabel('Longitude', 'fontsize', 14); ylabel('Latitude', 'fontsize', 14);
   
    %title(['slip-rate, cm/day, days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 12 );
    %title(['days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 10 );
    title([datestr1, ' - ',datestr2, ' ' ], 'fontsize', 10 );

    plot(coastline(:, 1), coastline(:, 2), 'k', 'Linewidth', 0.65)

    ylim(ylims)
    xlim(xlims)
    set(gca, 'fontsize', 10)
    if c > 1
        %set(h(c), 'YTickLabel', [], 'YTick', []);   
        set(h(c), 'YTickLabel', []);   
    else
        xlabel("Longtitude")
        ylabel("Latitude")
    end

    % Reduce the spacing by adjusting the subplot position
    pos = get(h(c), 'Position');
    pos(1) = pos(1) - (c-1) * 0.01; % Shift to the left; adjust the factor as needed
%     pos(2) = pos(2) - 0.01;
    pos(3) = pos(3) + 0.006;         % Increase width; adjust the factor as needed
    set(h(c), 'Position', pos);
end


c = 2*numplots;
for i = plot_list
    day1 = stationstruct(i).DOY;
    day2 = stationstruct(i+dur-1).DOY;
    year1=stationstruct(i).year;
    year2=stationstruct(i+dur).year;
    date1 = datetime(2011, 1, 1) + days(day1 - 1);  
    date2 = datetime(2011, 1, 1) + days(day2 - 1); 
    datestr1 = datestr(date1, 'mmm dd');
    datestr2 = datestr(date2, 'mmm dd');
 
    epoch = i;
    %slip rate
    slip_rate = mean(rate_b(:,epoch:(epoch+dur-1)), 2)*100/365;

    c=c+1;
   
    h(c) = subplot(1, 3*ceil(numplots), c);
    hold on

    %trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip_rate, 'EdgeAlpha', 0.225);

    validVertices = nd_ll(:,3) <= -20 & nd_ll(:,3) >= -60;
    validElements = all(validVertices(el), 2);
    filteredel = el(validElements, :);
    filteredsliprate = slip_rate(validElements, :);
    trisurf(filteredel,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),filteredsliprate, 'EdgeAlpha', 0.15);

    %caxis([-.25, .25])
    caxis([0, .25]);
    view(2);
    daspect([1/cosd(abs(origin(1))),1,110]);

    hold on;
    tremor_here = tremor(tremor(:,8) >= day1 & tremor(:,8) <= day2 & tremor(:,1) < 48.5, :);
    %plot(tremor_here(:,2), tremor_here(:,1), 'k.', 'MarkerSize', 4.4);
    
    %plot3(tremor_here(:,2), tremor_here(:,1), tremor_here(:,3), '.', 'MarkerSize', 0.01, 'Color', [1, 0.7, 0.85]);
    x = tremor_here(:, 2);
    y = tremor_here(:, 1);
    z = tremor_here(:, 3);
    scatter3(x, y, z, 2, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8 0 0.8], 'MarkerFaceAlpha', 0.2);

    %xlabel('Longitude', 'fontsize', 14); ylabel('Latitude', 'fontsize', 14);
   
    %title(['slip-rate, cm/day, days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 12 );
    %title(['days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 10 );
    title([datestr1, ' - ',datestr2, ' ' ], 'fontsize', 10 );

    plot(coastline(:, 1), coastline(:, 2), 'k', 'Linewidth', 0.65)

    ylim(ylims)
    xlim(xlims)
    set(gca, 'fontsize', 10)
    if c > 1
        %set(h(c), 'YTickLabel', [], 'YTick', []);   
        set(h(c), 'YTickLabel', []);   
    else
        xlabel("Longtitude")
        ylabel("Latitude")
    end
    
    % Reduce the spacing by adjusting the subplot position
    pos = get(h(c), 'Position');
    pos(1) = pos(1) - (c-1) * 0.01; % Shift to the left; adjust the factor as needed
%     pos(2) = pos(2) - 0.01;
    pos(3) = pos(3) + 0.006;         % Increase width; adjust the factor as needed
    set(h(c), 'Position', pos);
end

annotation('textbox', [0.15 0.66 0.22 0.02], 'String', 'Slip rates from FMOU', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none', ...
          "FontSize", 16);

annotation('textbox', [0.368 0.66 0.22 0.02], 'String', 'Slip rates from NIF', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none', ...
          "FontSize", 16);

annotation('textbox', [0.605 0.66 0.22 0.02], 'String', 'Slip rates from modified NIF', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none', ...
          "FontSize", 16);

 

% numColors = 512; 
% c = linspace(0, 1, numColors/4)'.^2.5; 
% whiteToBlue = [1-c, 1-c, ones(numColors/4, 1)];
% blueToGreen = [zeros(numColors/4, 1), linspace(0, 1, numColors/4)', ones(numColors/4, 1)];
% greenToYellow = [linspace(0, 1, numColors/4)', ones(numColors/4, 1), linspace(1, 0, numColors/4)'];
% yellowToRed = [ones(numColors/4, 1), linspace(1, 0, numColors/4)', 0.1+zeros(numColors/4, 1)];
% customCMap = [whiteToBlue; blueToGreen; greenToYellow; yellowToRed];

numColors = 900; 
c = linspace(0, 1, numColors/3)'.^2; 
whiteToBlue = [1-c, 1-c, ones(numColors/3, 1)];
blueToGreen = [zeros(numColors/3, 1), linspace(0, 1, numColors/3)', ones(numColors/3, 1)];
greenToYellow = [linspace(0, 1, numColors/3)', ones(numColors/3, 1), linspace(1, 0, numColors/3)'];
customCMap = [whiteToBlue; blueToGreen; greenToYellow];

% Apply the colormap
colormap(customCMap);

% Add a colorbar
cb = colorbar('Position',  [0.662 0.358 0.168 0.02], 'Orientation', 'horizontal');
% Set the colorbar's label
xlabel(cb, 'cm/day');
%set(gcf, 'Renderer', 'painters')
set(gcf, 'Position', [90, 90, 1250, 800]); 
%set(gcf, 'PaperOrientation', 'landscape');
%set(gcf, 'PaperPositionMode', 'auto');
%sgt = sgtitle('Slip rates from NIF', "FontSize", 15,'FontWeight', 'bold');

% print(gcf, '-dpng', '-r1200', '/Users/liuxubo/Desktop/real_data_FMOU_NIF_slip_rate_5days.png');

