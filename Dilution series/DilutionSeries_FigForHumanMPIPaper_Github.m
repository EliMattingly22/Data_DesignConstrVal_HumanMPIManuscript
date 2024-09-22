clear; clc; close all;
plot_on = 0;
set(0,'DefaultFigureWindowStyle','docked');

currdir = pwd;
addpath(strcat(pwd(),'\Data'))
% addpath(strcat(pwd(),'\Called_Functions'))
fflIter = 15;

%%
load('tSNRCompilation.mat')


FitSNR1_Signal0 = 0;
FFL1_IRad0 = 1;
UseEveryOther=1;
DetectionLimitSNR = 1;

conv_kernel_size_mm=6;
FOV_mm = 182;%mm


ImgSize = size(tSNRCompilation(1).fflImgsConcat(:,:,1));

res = FOV_mm/(ImgSize(1)-1);



conv_kernel_size_Pix_FWHM=conv_kernel_size_mm/res/(2*sqrt(2*log(2)));

stdcm_FFLImSpace =6;  %%mm

H = fspecial('gaussian', ImgSize,conv_kernel_size_Pix_FWHM);


Sig_irVec = [];
Sig_VecAll = [];
Sig_FFLVec = [];
SNR0_IRVec = [];
SNR0_IRVec_std = [];
SNR0_IRVec_Corner = [];
SNR0_IRVec_std_Corner = [];
SNR0_FFLVec = [];
SNR0_FFLVec_std = [];
SNR0_FFLVec_Corner = [];
SNR0_FFLVec_std_Corner = [];
RegressionMat_Sig = [0,0];
RegressionMat_SNR = [0,0];
massVec = [];
tSNR_IRVec = [];
load('tSNRAnalysis.mat')
load('tSNRCompilation.mat')



MassArray = [[tSNRCompilation.SPION_Mass_vec]', [1:length(tSNRCompilation)]'];

MassArraySort =sortrows(MassArray);





[~,BiggestSample] = max([tSNRCompilation.SPION_Mass_vec]);
if FFL1_IRad0
    [~,MaxSNRIndex] = max(tSNRAnalysis(BiggestSample).fflSNR0,[],'all');
    [row,col] = ind2sub(size(tSNRAnalysis(BiggestSample).fflSNR0),MaxSNRIndex);
else

    [~,MaxSNRIndex] = max(tSNRAnalysis(BiggestSample).irSNR0,[],'all');
    [row,col] = ind2sub(size(tSNRAnalysis(BiggestSample).irSNR0),MaxSNRIndex);
end
ROI_PlusMinus= 0;
ROI = ([row-ROI_PlusMinus ...
    row+ROI_PlusMinus ...
    col-ROI_PlusMinus ...
    col+ROI_PlusMinus]);

row2 = row+15;
col2 = col+15;
ROI2 = ([row2-ROI_PlusMinus ...
    row2+ROI_PlusMinus ...
    col2-ROI_PlusMinus ...
    col2+ROI_PlusMinus]);




H2 = fspecial('gaussian', ImgSize,conv_kernel_size_Pix_FWHM);

ImageSpace =zeros(ImgSize);
r = (ImgSize(1)-1)/2;
for j=1:ImgSize(1)
    for k=1:ImgSize(1)
        d = sqrt((j-0.5-(ImgSize(1)/2))^2 + (k-0.5-(ImgSize(1)/2))^2);
        if d<r
            ImageSpace(j,k) = 1;
        end
    end
end
ImageSpace = imfilter(ImageSpace,H2,'replicate');
% ImageSpace = ones(size(ImageSpace));
% tSNRAnalysis.irImgStack
%%
for kk=1:length(round(MassArraySort(2:end,2)))
    i = round(MassArraySort(1+kk,2));
    if FFL1_IRad0
        Sig = mean(tSNRAnalysis(i).fflImg_mean(ROI(1):ROI(2),ROI(3):ROI(4),:),'all');
        Sig_FFLVec = [Sig_FFLVec;Sig];
        SigAll = squeeze(mean(mean(tSNRAnalysis(i).fflImgStack(ROI(1):ROI(2),ROI(3):ROI(4),:),1),2));
        Sig_VecAll = [Sig_VecAll;SigAll];

        SNR0 = mean(tSNRAnalysis(i).fflSNR0(ROI(1):ROI(2),ROI(3):ROI(4),:),'all');
        SNR0Vec = squeeze(mean(mean(tSNRAnalysis(i).fflSNR0_Stack(ROI(1):ROI(2),ROI(3):ROI(4),:),1),2));
        SNR0_std = std(squeeze(tSNRAnalysis(i).fflSNR0_Stack(row,col,:)));
        SNR0_FFLVec = [SNR0_FFLVec;SNR0];
        SNR0_FFLVec_std = [SNR0_FFLVec_std;SNR0_std];

        SNR0_Corner = mean(tSNRAnalysis(i).fflSNR0(ROI2(1):ROI2(2),ROI2(3):ROI2(4),:),'all');
        SNR0_std_Corner = std(squeeze(tSNRAnalysis(i).fflSNR0_Stack(row2,col2,:)));
        SNR0_FFLVec_Corner = [SNR0_FFLVec_Corner;SNR0_Corner];
        SNR0_FFLVec_std_Corner = [SNR0_FFLVec_std_Corner;SNR0_std_Corner];
    else
        Sig = mean(tSNRAnalysis(i).irImg_mean(ROI(1):ROI(2),ROI(3):ROI(4),:),'all');
        Sig_irVec = [Sig_irVec;Sig];
        SigAll = squeeze(mean(mean(tSNRAnalysis(i).irImgStack(ROI(1):ROI(2),ROI(3):ROI(4),:),1),2));
        Sig_VecAll = [Sig_VecAll;SigAll];

        SNR0 = mean(tSNRAnalysis(i).irSNR0(ROI(1):ROI(2),ROI(3):ROI(4),:),'all');
        SNR0Vec = squeeze(mean(mean(tSNRAnalysis(i).irSNR0_Stack(ROI(1):ROI(2),ROI(3):ROI(4),:),1),2));
        SNR0_std = std(squeeze(tSNRAnalysis(i).irSNR0_Stack(row,col,:)));
        SNR0_IRVec = [SNR0_IRVec;SNR0];
        SNR0_IRVec_std = [SNR0_IRVec_std;SNR0_std];

        SNR0_Corner = mean(tSNRAnalysis(i).irSNR0(ROI2(1):ROI2(2),ROI2(3):ROI2(4),:),'all');
        SNR0_std_Corner = std(squeeze(tSNRAnalysis(i).irSNR0_Stack(row2,col2,:)));
        SNR0_IRVec_Corner = [SNR0_IRVec_Corner;SNR0_Corner];
        SNR0_IRVec_std_Corner = [SNR0_IRVec_std_Corner;SNR0_std_Corner];
    end

figure(23)
hold on
Hg = histogram(SigAll,3);
Hg.FaceColor = 'r';


mass = MassArraySort(end-i+1,1);

    RegressionMat_SNR = [RegressionMat_SNR;SNR0Vec,1000*mass*ones(size(SNR0Vec))];
    RegressionMat_Sig = [RegressionMat_Sig;SigAll,1000*mass*ones(size(SigAll))];
    massVec = [massVec;mass];


    numRows = 3;
    figure(45)
    
    
    
        ShowImgNum = 1;
        %     subplot(2,ceil((size(tSNRAnalysis,2)+1)/2),size(tSNRAnalysis,2)-i+1)
        subplot(numRows,ceil((size(tSNRAnalysis,2)+1)/numRows),size(tSNRAnalysis,2)-kk+1)

        if FFL1_IRad0
            imagesc( tSNRAnalysis(i).fflImgStack(:,:,ShowImgNum).*ImageSpace)
            %         imagesc(tSNRAnalysis(i).fflSNR0.*ImageSpace)
        else
            imagesc(tSNRAnalysis(i).irImgStack(:,:,ShowImgNum).*ImageSpace)
            %         imagesc(tSNRAnalysis(i).irSNR0.*ImageSpace)
        end

        Axis_mm = linspace(0,FOV_mm,4);
        Axis_mm_all = linspace(0,FOV_mm,66);
        xticks = linspace(1, 66, numel(Axis_mm));
        xticks_2 = linspace(1, FOV_mm, numel(Axis_mm));

        %     set(gca, 'XTick', xticks, 'XTickLabel', round(Axis_mm))
        %     set(gca, 'YTick', xticks, 'YTickLabel', round(Axis_mm))
        set(gca, 'XTick', [], 'XTickLabel', round(Axis_mm))
        set(gca, 'YTick', [], 'YTickLabel', round(Axis_mm))
        set(gca,'FontSize',13)
        axis image
        colormap hot
        ca = caxis();
        caxis([0 ca(2)]);
        
            if (tSNRCompilation(i).SPION_Mass_vec(1)*1000)<1000
                title(['',num2str(round(mass*1000)),'ng'],'FontSize',14,'FontWeight','bold')
            elseif (tSNRCompilation(i).SPION_Mass_vec(1))<10
                title(['',num2str(round(mass,1)),'\mu g'],'FontSize',14,'FontWeight','bold')
            else
                title(['',num2str(round(mass)),'\mu g'],'FontSize',14,'FontWeight','bold')
            end
       
   
    if i==1
        figure(45)
        subplot(numRows,ceil((size(tSNRAnalysis,2)+1)/numRows),ceil((size(tSNRAnalysis,2)+1)/numRows)*numRows)
        %     imagesc( log10(abs(mean(tSNRAnalysis(i).irImgStack(:,:,1:10),3))))
        if FFL1_IRad0
            EmptyImg = ((mean(tSNRCompilation(end).fflImgsConcat(:,:,ShowImgNum),3)));
        else
            EmptyImg = ((mean(tSNRCompilation(end).irImgsConcat(:,:,ShowImgNum),3)));
        end
        EmptyImg = imfilter(EmptyImg,H,'replicate').*ImageSpace;

        imagesc(EmptyImg)
        Axis_mm = linspace(0,FOV_mm,4);
        Axis_mm_all = linspace(0,FOV_mm,66);
        xticks = linspace(1, 66, numel(Axis_mm));
        xticks_2 = linspace(1, FOV_mm, numel(Axis_mm));

        %         set(gca, 'XTick', xticks, 'XTickLabel', round(Axis_mm))
        %         set(gca, 'YTick', xticks, 'YTickLabel', round(Axis_mm))
        set(gca, 'XTick', [], 'XTickLabel', round(Axis_mm))
        set(gca, 'YTick', [], 'YTickLabel', round(Axis_mm))
        set(gca,'FontSize',13)
        ca = caxis();
        caxis([0 ca(2)]);
        axis image
        colormap hot
        title('Empty','FontSize',14,'FontWeight','bold')
    end

    %     colorbar
    %     caxis([-8 -2])





end


if FFL1_IRad0
    UsedSNRData = SNR0_FFLVec;
    usedErrorData = SNR0_FFLVec_std;
    UsedSigData = Sig_FFLVec;
    UsedNoise = tSNRAnalysis(1).ffl_sigma0;

else
    UsedSNRData = SNR0_IRVec;
    usedErrorData = SNR0_IRVec_std;
    UsedSigData = Sig_irVec;
    UsedNoise = tSNRAnalysis(1).ir_sigma0;
end

massVec_ng = MassArraySort(MassArraySort(:,1)>0,1)*1000;
if FitSNR1_Signal0
    BestFit = [massVec_ng(:)\UsedSNRData(:), 0];
    MinDetectionMass = DetectionLimitSNR/BestFit(1);
    [xData, yData] = prepareCurveData( RegressionMat_SNR(:,2), RegressionMat_SNR(:,1) );

else
    BestFit = [massVec_ng(:)\UsedSigData(:), 0];
    MinDetectionMass = DetectionLimitSNR*UsedNoise/BestFit(1);
    [xData, yData] = prepareCurveData( RegressionMat_Sig(:,2), RegressionMat_Sig(:,1) );

end



% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf 0];
opts.Upper = [Inf 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

BestFitTestMasses = [MinDetectionMass/5;MinDetectionMass;massVec_ng(:)];
BestFit_Vals = polyval([fitresult.p1 fitresult.p2],BestFitTestMasses);


figure(20),
hold on

%         errorbar(massVec_ng,UsedSNRData,usedErrorData,'bd','LineWidth',2)
if FitSNR1_Signal0
        plot(RegressionMat_SNR(:,2),RegressionMat_SNR(:,1),'b*','LineWidth',2)
    loglog(BestFitTestMasses,BestFit_Vals,'--','LineWidth',2,'Color',[0.2,0.2,.2])
    loglog(MinDetectionMass,polyval(BestFit,MinDetectionMass),'ro','LineWidth',2)
    ylabel('SNR','FontSize',14,'FontWeight','bold')
    legend({'Measured Data',['Best Fit, R^2=',num2str(gof.rsquare)],'Predicted Detection Limit'},'Location','northwest')


else
    plot(RegressionMat_Sig(:,2),RegressionMat_Sig(:,1),'b*','LineWidth',2)
    loglog(BestFitTestMasses,BestFit_Vals,'--','LineWidth',2,'Color',[0.2,0.2,.2])
    loglog(MinDetectionMass*5,polyval(BestFit,MinDetectionMass*5),'ro','LineWidth',2)
    
    loglog([min(BestFitTestMasses) 5*max(BestFitTestMasses)],[UsedNoise UsedNoise],'-','LineWidth',2,'Color',[0.2,0.5,.5])
    loglog([min(BestFitTestMasses) 5*max(BestFitTestMasses)],5*[UsedNoise UsedNoise],'--','LineWidth',2,'Color',[0.2,0.5,.5])
    
    ylabel('Signal [A.U]','FontSize',14,'FontWeight','bold')
    legend({'Measured Data',['Best Fit, R^2=',num2str(gof.rsquare,3)],'SNR=5 Detection Limit','1\sigma Noise Floor','5\sigma Noise Floor'},'Location','northwest')

end



grid on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('Mass Synomag 70nm [ng Fe]','FontSize',14,'FontWeight','bold')

if FFL1_IRad0 == 1
    title('FFL Recon')
else
    title('Inverse Radon Recon')
end



    xlim([100 2e5])
%     ylim([1,3e4])

set(gca,'FontSize',12)


if FFL1_IRad0
    NoiseCollection = imfilter(tSNRCompilation(end).fflImgsConcat,H2,'replicate');
else
    NoiseCollection = imfilter(tSNRCompilation(end).irImgsConcat,H2,'replicate');
end
NoiseCollection = NoiseCollection(20:40,20:40,:);
NoiseCollection = NoiseCollection(:);

% 
% figure(23),histogram(NoiseCollection,70)
% StrMeanNoise =  num2str(mean(NoiseCollection));
% StrStdNoise = num2str(std(NoiseCollection));
% title({'Histogram of signal of pixels in empty bore',['Mean =',StrMeanNoise,'Std = ',StrStdNoise]})
% % RSqr(RegressionMat,BestFit)



function [R2]=RSqr(RegMat,FitParam)
RegMatWithFit = [RegMat,RegMat(:,2).*FitParam];
SSE = sum((RegMatWithFit(:,1)-RegMatWithFit(:,3)).^2);
SST = sum((RegMatWithFit(:,1)).^2);
R2 = 1-SSE/SST;
end

function [X,Y] = prepareCurveData(X,Y)
X = X(:);
Y = Y(:);
end