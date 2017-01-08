%% EECS 442 - HW 04 - Q2 Blob Detection
%  ------------------------------------------------------------------------
%  Date: 11 / 21 / 2016
%  Author: Fan Bu
%  ------------------------------------------------------------------------
%  show_all_circles(I, cx, cy, rad, color, ln_wid)
%  ------------------------------------------------------------------------
%  Instructions:
%  Implement a blob detector based on the Laplacian operator and 
%  characteristic scale concept introduced by Lindeberg. Compute the radius
%  of each blob (in pixels) and report the histogram (dis- tribution)
%  of radii computed for the whole image.
% 
% Algorithm outline:\\
% 1. Build a Laplacian scale space, starting with some initial scale and going for $n$ iterations.\\
% 2. Generate a scale-normalized Laplacian of Gaussian filter at a given scale "$sigma$".\\
% 3. Filter image with the scale-normalized Laplacian. \\
% 4. Save square of Laplacian filter response for current level of scale space. \\
% 5. Determine threshold by multiply the scale with a factor $k$. \\
% 6. Perform non-maximum suppression in scale space. \\
% 7. Display resulting circles at their characteristic scales.\\

%% Initialization
clear all;
close all;
clc;
%% ------------------- Parameter selection -----------------------
t = 1:0.1:3;
sigma = exp(t);
sigma = sigma';
%levels = 10;
levels = length(sigma);
%% ------------------- Load data --------------------------------
img = imread('3.bmp');
img = rgb2gray(img);
double_img = double(img);

% Build a Laplacian scale space, starting with some initial scale and going for $n$ iterations.
[height width] = size(img);
scale_space = zeros(height,width,levels);

%% ------------------- Perform detection ------------------------

for i=1:levels
    
% Generate a scale-normalized Laplacian of Gaussian filter at a given scale "$sigma$".
filter_size = 2*ceil(3*sigma(i))+1;
kernel = fspecial('log', filter_size, sigma(i));
% Filter image with the scale-normalized Laplacian.
nLoG = sigma(i)^2 * kernel;
Filtered_img = imfilter(double_img, nLoG, 'same', 'replicate');
Filtered_img = Filtered_img .^ 2;
% Save square of Laplacian filter response for current level of scale space.
scale_space(:,:,i) = Filtered_img;
end
% choose the right thresold T, as T = k * scaleSpace_{max}.
k = 0.4;
max(max(max(scale_space)))
T = k * max(max(max(scale_space)))

% Perform non-maximum suppression in the 3D scale space
scale_space_filtered = zeros(height,width,levels);
for i=1:levels
    scale_space_filtered(:,:,i) = ordfilt2(scale_space(:,:,i),9,ones(3));
end
for i = 1:levels
    scale_space_filtered(:,:,i) =  max(scale_space_filtered(:,:,max(i-1,1):min(i+1,levels)),[],3);
end
scale_space_filtered  = scale_space_filtered .* (scale_space_filtered == scale_space);

% calculate perameters for show_all_circles
row = [];
col = [];
radius = [];
for i=1:levels
    [new_rows,new_cols] = find(scale_space_filtered(:,:,i)>=T);
    num = length(new_rows);
    radii = sigma(i)*sqrt(2);
    radii = repmat(radii,num,1);
    row = [row;new_rows];
    col = [col;new_cols];
    radius = [radius;radii];
end
show_all_circles(img,col,row,radius);
set(gcf,'units','normalized','position',[0,0,0.8,1]);
h_title=title({['Number of circles = ', num2str(size(col,1))]
    ['Safety factor k = ',num2str(k), ';  Threshold T = ',num2str(T)]});
set(h_title,'FontSize',18);
print(gcf,'-djpeg' ,['HW4_q2_k_',num2str(k),'.jpeg'],'-r400')