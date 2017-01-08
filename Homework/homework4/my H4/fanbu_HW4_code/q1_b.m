%% EECS 442 - HW 04 - Q1 Harris Corner Detection
%  ------------------------------------------------------------------------
%  Date: 11 / 21 / 2016
%  Author: Fan Bu
%  ------------------------------------------------------------------------
%  ellipse(ra,rb,ang,x0,y0,C,Nb)
%  ------------------------------------------------------------------------
%  Instructions(Latex code)
% Key steps of the algorithm:\\
% Step1, Compute magnitude of the x and y gradients at each pixel.\\
% Step2, Due to Harris, construct a Gaussian window M around each pixel.\\
% Step3, Compute $\lambda_s$ of second moment matrix $M$ \\
% Step4, Compute $R=det(M)-k(tr(M))^2$ \\
% Step5, Set proper threshold $T$ \\
% Step6, If $R> T$, a corner is detected; then, retain point of local maxima. \\
% Step7, Superpose corners on images.\\

%% Initialization
clear all; 
close all; 
clc;
%% ----------------------- Load data ----------------------------
img1 = imread('I1.png'); 
img2 = imread('I3.png');

img = img1;
w = 9; %choose window size

%% ------------ Part 1: compute M,R and ellipse parameter --------
% translate to gray image
gray_img = rgb2gray(img);
[height,width] = size(gray_img);   
% initialize corner position
is_edge = zeros(height, width); 
% initialize R score
R = zeros(height, width);

% applying sobel edge detector in the horizontal direction
fx = [-1 0 1;
    -1 0 1;
    -1 0 1];
Ix = filter2(fx,gray_img);
% applying sobel edge detector in the vertical direction
fy = [1 1 1;
    0 0 0;
    -1 -1 -1];
Iy = filter2(fy,gray_img); 
% compute terms for second moment matrix M
Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;

% perform gaussian filtering for eliminating noises
h= fspecial('gaussian',[w w],2); 
Ix2 = filter2(h,Ix2);
Iy2 = filter2(h,Iy2);
Ixy = filter2(h,Ixy);

% Initialize ellipsePerem
ellipsePerem = [];
Rmax = 0;
figure;
imshow(img(1:height/2,1:width/2,:),'border','tight','initialmagnification','fit');
hold on;
for i = 1:height
    for j = 1:width
        % compute second moment matrix M
        M = [
            Ix2(i,j), Ixy(i,j);
            Ixy(i,j) Iy2(i,j)
            ]; 
        % compute R(i,j) = det(M)-0.04*(trace(M))^2
        R(i,j) = det(M)-0.04*(trace(M))^2;
        if R(i,j) > Rmax
           Rmax = R(i,j);
        end;

    end;
end;
num_edge = 0;
Rmax
T=0.1*Rmax
% determine corner by comparing with threshold and all neighbor points
for i = 2:height-1
    for j = 2:width-1
        M = [
                 Ix2(i,j), Ixy(i,j);
                 Ixy(i,j) Iy2(i,j)
            ]; 
        % compute M's eigenvalue
        [vector,lambda] = eig(M);
        lambda = diag(lambda);
        % compute radius associate to the ellipse
        [lambda1,index] = max(lambda); 
        lambda2 = min(lambda); 
        a = lambda1^(-0.5);
        b = lambda2^(-0.5);
        vector = vector(:,index);
        angle = atan(vector(2)/vector(1));
        % form an ellipsePerem
        if mod(i,4) == 0 && mod(j,4) == 0
            ellipsePerem = [ellipsePerem;a b angle j i];
        end
        if R(i,j) > T ...
            && R(i,j) > R(i-1,j-1) ...
            && R(i,j) > R(i-1,j) ...
            && R(i,j) > R(i-1,j+1) ...
            && R(i,j) > R(i,j-1) ...
            && R(i,j) > R(i,j+1) ...
            && R(i,j) > R(i+1,j-1) ...
            && R(i,j) > R(i+1,j) ...
            && R(i,j) > R(i+1,j+1)
                
            is_edge(i,j) = 1;
            num_edge = num_edge+1;
        end;
    end;
end;
[edge_col, edge_row] = find(is_edge == 1);
num_edge
%% ========== Part 2: Plot ellipse and corner on image ====================

ellipse(ellipsePerem(:,1), ellipsePerem(:,2), ellipsePerem(:,3), ...
    ellipsePerem(:,4),ellipsePerem(:,5));
print(gcf,'-djpeg' ,['HW4_q1_b_w_',num2str(w),'_I1.jpeg'],'-r400')