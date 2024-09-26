clear all 
close all

load('TxHeatingData40A.mat')
 figure,
    plot(TimeSec0/60,Txmag_ij,'r-','LineWidth',3)
     set(gca,'FontSize',12)
    xlabel('Experiment time [Minutes]','FontSize',16,'FontWeight','Bold')
    ylabel({'Drive coil current [A_{pk}]'},'FontSize',16,'FontWeight','Bold')

    yl = ylim;
    ylim(yl)
    yyaxis right
    p2 = plot(TimeSec0/60,Txmag_ij,'r-','LineWidth',3);
    ax = gca;
    ylim(yl);
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
labelVals = get(ax,'YTick');
for kk =1:length(labelVals)
LabelCellVec(kk) = {num2str(100*(labelVals(kk)/mean(Txmag_ij)-1),2)};
end
set(gca,'yticklabels',LabelCellVec)
ylabel('Deviation from mean [%]','FontSize',16,'FontWeight','Bold')
    set(gcf,'Position',[4.2703e+03 133 936.6667 420])
    set(gca,'FontSize',14)
