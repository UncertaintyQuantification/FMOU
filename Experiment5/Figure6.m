load ../real_data/MatLab_code_Modified_NIF/example/Cascadia_xubo_setup4_f1b1W1_.mat;
load ../real_data/MatLab_code_Modified_NIF/example/GPS_info.mat;
load ../real_data/MatLab_code_Modified_NIF/example/coastline.dat;


x_center = -123.2;
y_center = 46.0;
semi_minor = 0.7;
semi_major = semi_minor/3;
v= 8;

FMOU_ellipse_slip = readmatrix("saved_data/FMOU_ellipse_slip_large_mesh.csv");
true_ellipse_slip = readmatrix("saved_data/true_ellipse_slip_large_mesh.csv");
NIF_ellipse_slip = readmatrix("saved_data/NIF_ellipse_slip_large_mesh.csv");


beginslip = 3; 
endslip = 70; 
skip = 12;
dur = 5;
xlims = [-125.6 -121.8];
ylims = [43.8 48.9];

plot_list = [3, 36, 69];
figure
c = 0;
%numplots=ceil((endslip-beginslip)/skip);
numplots=length(plot_list);
h = zeros(1, 3*numplots);

for i = plot_list
    day1 = stationstruct(i).DOY;
    day2 = stationstruct(i+dur-1).DOY;
    year1=stationstruct(i).year;
    year2=stationstruct(i+dur).year;
    date1 = datetime(2011, 1, 1) + days(day1 - 1);  
    date2 = datetime(2011, 1, 1) + days(day2 - 1); 
    datestr1 = datestr(date1, 'mmm dd');
    datestr2 = datestr(date2, 'mmm dd');
    % datestr1 = datestr(date1, 'mmmm dd'); % full name of each month
    % datestr2 = datestr(date2, 'mmmm dd');

    %slip rate
    slip = mean(FMOU_ellipse_slip(:,i:(i+dur-1)), 2);
    %integrated slip rate
    %slip = sum(rate_b(:,1:i), 2)/800;

    c=c+1;

    h(c) = subplot(1, 3*ceil(numplots), c);
    hold on;

    %trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip, 'EdgeAlpha', 0.225);
    
    validVertices = nd_ll(:,3) <= -20 & nd_ll(:,3) >= -60;
    validElements = all(validVertices(el), 2);
    filteredel = el(validElements, :);
    filteredsliprate = slip(validElements, :);
    trisurf(filteredel,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),filteredsliprate, 'EdgeAlpha', 0.225);

%     tremor_here = tremor(tremor(:,8) >= day1 & tremor(:,8) <= day2, :);
%     % plot(tremor_here(:,2), tremor_here(:,1), 'k.', 'MarkerSize', 4.4);
%     plot3(tremor_here(:,2), tremor_here(:,1), tremor_here(:,3), 'k.', 'MarkerSize', 3);
     
    caxis([0, 3.0]); 
    view(2);
    daspect([1/cosd(abs(origin(1))),1,110]);

    hold on;

    semi_major = semi_minor/3 + v*(year2 - stationstruct(1).year);
    position = [x_center-semi_minor, y_center-semi_major, 2 * semi_minor, 2 * semi_major];
    rectangle('Position', position, 'Curvature', [1, 1], 'EdgeColor', 'red','Linewidth', 2.25);

  
    %title(['slip-rate, cm/day, days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 12 );
    %title(['days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 10 );
    title([datestr1, ' - ',datestr2, ' ' ], 'fontsize', 9 );

    plot(coastline(:, 1), coastline(:, 2), 'k', 'Linewidth', 0.8)


    ylim(ylims)
    xlim(xlims)
    set(gca, 'fontsize', 10);    
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
    % datestr1 = datestr(date1, 'mmmm dd'); % full name of each month
    % datestr2 = datestr(date2, 'mmmm dd');

    %slip rate
    slip = mean(NIF_ellipse_slip(:,i:(i+dur-1)), 2);
    %integrated slip rate
    %slip = sum(rate_b(:,1:i), 2)/800;

    c=c+1;

    h(c) = subplot(1, 3*ceil(numplots), c);
    hold on;

    %trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip, 'EdgeAlpha', 0.225);
    
    validVertices = nd_ll(:,3) <= -20 & nd_ll(:,3) >= -60;
    validElements = all(validVertices(el), 2);
    filteredel = el(validElements, :);
    filteredsliprate = slip(validElements, :);
    trisurf(filteredel,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),filteredsliprate, 'EdgeAlpha', 0.225);

%     tremor_here = tremor(tremor(:,8) >= day1 & tremor(:,8) <= day2, :);
%     % plot(tremor_here(:,2), tremor_here(:,1), 'k.', 'MarkerSize', 4.4);
%     plot3(tremor_here(:,2), tremor_here(:,1), tremor_here(:,3), 'k.', 'MarkerSize', 3);
     
    caxis([0, 3.0]); 
    view(2);
    daspect([1/cosd(abs(origin(1))),1,110]);

    hold on;

    semi_major = semi_minor/3 + v*(year2 - stationstruct(1).year);
    position = [x_center-semi_minor, y_center-semi_major, 2 * semi_minor, 2 * semi_major];
    rectangle('Position', position, 'Curvature', [1, 1], 'EdgeColor', 'red','Linewidth', 2.25);

  
    %title(['slip-rate, cm/day, days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 12 );
    %title(['days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 10 );
    title([datestr1, ' - ',datestr2, ' ' ], 'fontsize', 9 );

    plot(coastline(:, 1), coastline(:, 2), 'k', 'Linewidth', 0.8)


    ylim(ylims)
    xlim(xlims)
    set(gca, 'fontsize', 10);    
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
    %pos(2) = pos(2) + 0.012;
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
    % datestr1 = datestr(date1, 'mmmm dd'); % full name of each month
    % datestr2 = datestr(date2, 'mmmm dd');

    %slip rate
    slip = mean(true_ellipse_slip(:,i:(i+dur-1)), 2);
    %integrated slip rate
    %slip = sum(rate_b(:,1:i), 2)/800;

    c=c+1;

    h(c) = subplot(1, 3*ceil(numplots), c);
    hold on;

    %trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip, 'EdgeAlpha', 0.225);
    
    validVertices = nd_ll(:,3) <= -20 & nd_ll(:,3) >= -60;
    validElements = all(validVertices(el), 2);
    filteredel = el(validElements, :);
    filteredsliprate = slip(validElements, :);
    trisurf(filteredel,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),filteredsliprate, 'EdgeAlpha', 0.225);

%     tremor_here = tremor(tremor(:,8) >= day1 & tremor(:,8) <= day2, :);
%     % plot(tremor_here(:,2), tremor_here(:,1), 'k.', 'MarkerSize', 4.4);
%     plot3(tremor_here(:,2), tremor_here(:,1), tremor_here(:,3), 'k.', 'MarkerSize', 3);
     
    caxis([0, 3.0]); 
    view(2);
    daspect([1/cosd(abs(origin(1))),1,110]);

    hold on;

    semi_major = semi_minor/3 + v*(year2/2 + year1/2 - stationstruct(1).year);
    position = [x_center-semi_minor, y_center-semi_major, 2 * semi_minor, 2 * semi_major];
    rectangle('Position', position, 'Curvature', [1, 1], 'EdgeColor', 'red','Linewidth', 2.6);

  
    %title(['slip-rate, cm/day, days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 12 );
    %title(['days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 10 );
    title([datestr1, ' - ',datestr2, ' ' ], 'fontsize', 9 );

    plot(coastline(:, 1), coastline(:, 2), 'k', 'Linewidth', 0.8)

    ylim(ylims)
    xlim(xlims)
    set(gca, 'fontsize', 10);
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
    %pos(2) = pos(2) + 0.032;
    pos(3) = pos(3) + 0.006;         % Increase width; adjust the factor as needed
    set(h(c), 'Position', pos);
end


annotation('textbox', [0.15 0.66 0.22 0.02], 'String', 'Slips from FMOU', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none', ...
          "FontSize", 16);

annotation('textbox', [0.368 0.66 0.22 0.02], 'String', 'Slips from NIF', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none', ...
          "FontSize", 16);

annotation('textbox', [0.605 0.66 0.22 0.02], 'String', 'True slips', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none', ...
          "FontSize", 16);

numColors = 512; 
c = linspace(0, 1, numColors/4)'.^2.5; 
whiteToBlue = [1-c, 1-c, ones(numColors/4, 1)];
blueToGreen = [zeros(numColors/4, 1), linspace(0, 1, numColors/4)', ones(numColors/4, 1)];
greenToYellow = [linspace(0, 1, numColors/4)', ones(numColors/4, 1), linspace(1, 0, numColors/4)'];
yellowToRed = [ones(numColors/4, 1), linspace(1, 0, numColors/4)', zeros(numColors/4, 1)];

customCMap = [whiteToBlue; blueToGreen; greenToYellow; yellowToRed];

% Apply the colormap
colormap(customCMap);

% Add a colorbar
cb = colorbar('Position',  [0.662 0.358 0.168 0.02], 'Orientation', 'horizontal');
% Set the colorbar's label
xlabel(cb, 'cm');
%set(gcf, 'Renderer', 'painters')
set(gcf, 'Position', [90, 90, 1250, 800]);
set(gcf, 'PaperOrientation', 'landscape');
%set(gcf, 'PaperPositione', 'auto');

%sgt = sgtitle('Elliptical Slips region, slip rates averaged over 3 days', "FontSize", 20,'FontWeight', 'bold');

%print(gcf, '-dpng', '-r600', 'FMOU_NIF_true_ellipse_slip.png'); 
