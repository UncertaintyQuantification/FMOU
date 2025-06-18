function slip_hist(beginslip, endslip, skip, dur, xlims, ylims, NIFsavefile, coastline)	
%slip_hist(beginslip, endslip, skip, dur, xlims, ylims, NIFsavefile, coastline)	
% plot the slip history and total slip output from the NIF
% inputs:
% beginslip: starting epoch index
% endslip: ending epoch index
% skip:epochs between panels of slip-rate snapshots
% dur: averaging interval for panels of slip-rate snapshots
% xlims: longitude limits of plot, 1 by 2 vector
% ylims: latitude limits of plot, 1 by 2 vextor
% NIFsavefile: string name of NIF output file
% coastline: array that contains coastline for plots, first col longitude points, second col latitude points

figure

c=0;

load GreyNegMapbig
load(NIFsavefile)

numplots=ceil((endslip-beginslip)/skip);

for i =beginslip:skip:endslip;

    day1 = stationstruct(i).DOY;
    day2 = stationstruct(i+dur-1).DOY;
    year1=stationstruct(i).year;
    year2=stationstruct(i+dur).year;

    %slip rate
    slip = mean(rate_b(:,i:(i+dur-1)), 2)*100/365;
    %integrated slip rate
    %slip = sum(rate_b(:,1:i), 2)/800;

    c=c+1;

    subplot(2, ceil(numplots/2), c);
    disp(c);
    hold on

    trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip, 'EdgeAlpha', 0.1);
    caxis([-.25, .25])
    view(2);
    daspect([1/cosd(abs(origin(1))),1,110]);

    %xlabel('Longitude', 'fontsize', 14); ylabel('Latitude', 'fontsize', 14);
    title(['slip-rate, cm/day, days ', num2str(day1), '-',num2str(day2), ' ' ], 'fontsize', 12 );

    plot(coastline(:, 1), coastline(:, 2), 'k', 'Linewidth', 2)


    ylim(ylims)
    xlim(xlims)
    set(gca, 'fontsize', 16)
end

set(gcf, 'Colormap', GreyNegMapbig);
colorbar
set(gcf, 'Renderer', 'painters')
  

    
    %%
    
figure

load GreyNegMapbig
slip = slip_b(:,end-5)-slip_b(:, 5);


hold on

trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip*100, 'EdgeAlpha', 0.2);

caxis([-2.5, 2.5])
view(2);
daspect([1/cosd(abs(origin(1))),1,110]);

set(gcf, 'Colormap', GreyNegMapbig);
    
plot(coastline(:, 1), coastline(:, 2), 'k', 'Linewidth', 2)

   
ylim(ylims)
xlim(xlims)
set(gca, 'fontsize', 16)
set(gcf, 'Renderer', 'painters')
colorbar
title('Total Slip, cm', 'fontsize', 14)

    
   
