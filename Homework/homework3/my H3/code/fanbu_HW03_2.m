function main()
clear all;
close all;
clc;
%% Load in images
image1 = imread('1.JPG');
image2 = imread('2.JPG');
figure, imagesc(image1); colormap(gray);
title('Click a point on this image'); axis image;
%% Select corner - make sure to click in order corresponding to A
%% If you want to use my data, just comment this subblock
% num_of_points=12;
% a = [];
% for i = 1:num_of_points
%     [x y] = ginput(1);
%     a = [a [x;y]]
%     image1(round(y),round(x)) = 255;
%     imagesc(image1); colormap(gray);
%     title('Click a point on this image');
%     axis image;
% end
% figure, imagesc(image2); colormap(gray);
% title('Click a point on this image');
% axis image;
% for i = 1:num_of_points
%     [x y] = ginput(1);
%     a = [a [x;y]]
%     image1(round(y),round(x)) = 255;
%     imagesc(image2); colormap(gray);
%     title('Click a point on this image');
%     axis image;
% end
% save ground_truth
%% Load my data
load ground_truth.mat
%% 3D coordinates
%Build the reference A matrix
A = [0 0 0]';
for i=1:5
    A = [A A(:,end)];
    A(1,end)=A(1,end)+35;
end
A = [A [0 35 0]'];
for i=1:5
    A = [A A(:,end)];
    A(1,end)=A(1,end)+35;
end
A = [A [0 0 9.5]'];
for i=1:5
    A = [A A(:,end)];
    A(1,end)=A(1,end)+35;
end
A = [A [0 35 9.5]'];
for i=1:5
    A = [A A(:,end)];
    A(1,end)=A(1,end)+35;
end
A = [A;ones(1,size(A,2)) ]
A_annot = [a; ones(1, size(a,2))]
%Solve using least squares
P = Calc_P(A', a')
%RMS Calculation
A_bar = P * A
RMS_error = sqrt(mean(sum((A_bar - A_annot).^2,1)))
% RMS = 1.3175
%--------------------------------------------------------
%returns the H matrix given point (3D) and pointp(image) (transformed coordinates)
% point is N X 3, pointp is N X 2

end

function P = Calc_P (standard_point, annot_point)
A = [];
b = [];
%Build matrix A and b
for i = 1:size(standard_point,1)
    x_now = standard_point(i,1);
    x_anno_now = annot_point(i,1);
    y_now = standard_point(i,2);
    y_anno_now = annot_point(i,2);
    z_now =standard_point(i, 3);
    %calculates A matrix and B vector
    A = [A;
        x_now y_now z_now 1 0 0 0 0;
        0 0 0 0 x_now y_now z_now 1];
    b = [b; x_anno_now; y_anno_now];
end
%least squares solution
x = (A'*A)^-1 * A'*b;
P = [x(1) x(2) x(3) x(4);
    x(5) x(6) x(7) x(8);
    0 0 0 1];
end