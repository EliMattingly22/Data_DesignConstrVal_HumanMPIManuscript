%% Importing all h5 Files and sort Image data
clear all
close all
ptsBG = 50;
set(0,'DefaultFigureWindowStyle','normal')
load('CollectedResData_2VShift.mat')


Nfiles = size(FFLimgCompile,3);

UseIter = 15;
ContrastVal_IR = zeros(Nfiles,1);
DriveAmp = zeros(Nfiles,1);
PkRatio_FFL = zeros(Nfiles,1);
figure(99), clf
for i = 1 : Nfiles
    
    
    
    kernelFWHM_mm = 1;
    kernelFWHM_pix = kernelFWHM_mm/FOV_img*ROPPs;
    kernel_std = round(kernelFWHM_pix/2.355);%FWHM to Sigma

    
    IRimg = imgaussfilt(IRimgCompile(:,:,i),kernel_std);

     
    FFLimg = abs(imgaussfilt(FFLimgCompile(:,:,i),kernel_std));

    
    subplot(3,Nfiles,i)
    imagesc(IRimg)
    axis image
    % colormap inferno
    hold on
    plot([1 ROPPs],[ROPPs/2 ROPPs/2],'g-','LineWidth',2)
    if i==1
        ylabel({'IRadon Recon'},'interpreter','none','FontSize',14,'FontWeight','bold')
    end
    title(['Separation = ',num2str(dist(i)),' mm'],'FontSize',14,'FontWeight','bold')
    set(gca,'XTickLabel',{},'YTickLabel',{})


    subplot(3,Nfiles,2*Nfiles+i), hold on
    XVals = linspace(-0.5,0.5,ROPPs)*FOV_img;
    YVals = IRimg(ROPPs/2,:);
    YVals = YVals/max(YVals);
    plot(XVals,YVals,'g-','LineWidth',2)

    [PkVals,PkLocs] = findpeaks(YVals,"NPeaks",2,"MinPeakDistance",5,'MinPeakProminence',max(YVals)/15,'SortStr','descend');
    TroughVal = min(YVals(min(PkLocs):max(PkLocs)));
    ContrastVal_IR(i) = 1-TroughVal/mean(PkVals);
    xlim([min(XVals) max(XVals)])

    subplot(3,Nfiles,Nfiles+i)
    imagesc(FFLimg)
    axis image
    % colormap inferno
    hold on
    plot([1 ROPPs],[ROPPs/2 ROPPs/2],'c:','LineWidth',2)
    if i==1
        ylabel({'Iterative Recon'},'interpreter','none','FontSize',14,'FontWeight','bold')
    end

    set(gca,'XTickLabel',{},'YTickLabel',{})
  subplot(3,Nfiles,2*Nfiles+i),hold on
    XVals = linspace(-0.5,0.5,ROPPs)*FOV_img;
    YVals = FFLimg(ROPPs/2,:);
    YVals = YVals/max(YVals);
    plot(XVals,YVals,'c:','LineWidth',2)
    if i==1
        ylabel({'Normalized';'Signal [A.U]'},'FontSize',14,'FontWeight','bold')
    end
    xlabel('Dist. along line [mm]','FontSize',14,'FontWeight','bold')
    set(gca,'FontSize',12)
    xlim([min(XVals) max(XVals)])
    [PkVals,PkLocs] = findpeaks(YVals,"NPeaks",2,"MinPeakDistance",5,'MinPeakProminence',max(YVals)/15,'SortStr','descend');
    TroughVal = min(YVals(min(PkLocs):max(PkLocs)));
    ContrastVal(i) = 1-TroughVal/mean(PkVals);

end

figure,
plot(dist,ContrastVal_IR,'g*','LineWidth',2)
hold on
plot(dist,ContrastVal,'cd','LineWidth',2)
xlabel('Separation [mm]','FontSize',16,'FontWeight','bold')
ylabel('Signal Contrast [UL]','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',14)
legend('Inv. Radon','Iterative Recon','Location','northwest')