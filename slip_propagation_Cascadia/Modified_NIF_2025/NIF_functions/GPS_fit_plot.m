function GPS_fit_plot(plotlist, NIFsavefile)
%GPS_fit_plot(plotlist, NIFsavefile)
% plot GPS data and fits output from the NIF
% east data in blue, north data offset in red
% inputs:
% plotlist: a 4 by n character array of station names
% NIFsavefile: string name of NIF output file

load(NIFsavefile)

days=[stationstruct(:).DOY]-stationstruct(1).DOY+1;
days=[stationstruct(:).year];

numplots=ceil(length(plotlist));

figure
c=0;
for j=1:length(plotlist)
for i=1:length(ID)
    
    if(GetIndex(plotlist(j, :), ID(i, 1:4))>0)
        
        subplot(ceil(numplots/2), 2, j);
    title(ID(i,1:4), 'fontsize', 14)
    hold on
    plot(days,E(i, 1:end)*100-TT(1, :)*100, '.', 'MarkerEdgeColor', [.15 .3 1])
    plot(days, Ep(i, :)*100, 'k', 'LineWidth', 2);
    
    plot(days,N(i, 1:end)*100-1-TT(2, :)*100, 'r.')
    plot(days, Np(i, :)*100-1, 'k', 'LineWidth', 2);
    
    minN=min(N(i, :))*100-1;
    
    %uncomment these lines to add verticals to plot
    %plot(days,U(i, 1:end)*100+2*minN, 'g.')
    %plot(days, Up(i, :)*100+2*minN, 'k', 'LineWidth', 2);
    
    axis tight
    set(gca, 'fontsize', 14)
    xlabel('Year', 'fontsize', 14)
    ylabel('Position, cm', 'fontsize', 14)
    
    end

   
   
end
end

