
clear; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');
figs_on = 0;

load('2023_10_27_14_04_32_Big-G-phantom_51Apk_5Vshift_pp_FFL_irad.mat');
% size is 66 x 66 x 20 iters 

iter = 25;
% img = FFL_recon.image000.images(:,:,iter)+FFL_recon.image002.images(:,:,iter)+FFL_recon.image004.images(:,:,iter);
img = iradon_recon.image000.IR+iradon_recon.image002.IR+iradon_recon.image004.IR;
img = imgaussfilt(flipud(img),2);

ImgSize = size(img);
H2 = fspecial('gaussian', ImgSize,5);

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

img = img.* ImageSpace;

figure(68),
% subplot(1,2,1)
imagesc(img); colormap jet; axis image;

phantom_dist = (132+140)/2; % [mm]
pt1 = [69, 13];

pt2 = [71, 112];

hold on, plot(pt1(1),pt1(2),'*','Color','blue')
hold on, plot(pt2(1),pt2(2),'*','Color','blue')
colormap hot

c = caxis;
caxis([0 c(2)]);


ydist_pixels = pt2(2)-pt1(2);
xdist_pixels = pt2(1)-pt1(1);
dist_pixels = sqrt(xdist_pixels^2 + ydist_pixels^2);

mm_per_pixel = phantom_dist/dist_pixels;

FOV_mm = mm_per_pixel*size(img,1); % in mm!! 

FOV_m = FOV_mm*1e-3;

title(['FOV =',num2str(FOV_mm,3),' mm'])
  Axis_mm = linspace(0,FOV_mm,4);
        Axis_mm_all = linspace(0,FOV_mm,66);
        xticks = linspace(1, 66, numel(Axis_mm));
        xticks_2 = linspace(1, FOV_mm, numel(Axis_mm));

        %     set(gca, 'XTick', xticks, 'XTickLabel', round(Axis_mm))
        %     set(gca, 'YTick', xticks, 'YTickLabel', round(Axis_mm))
        set(gca, 'XTick', [], 'YTick', [])
        set(gca,'FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
