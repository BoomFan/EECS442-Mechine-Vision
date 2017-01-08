function main
clear all;
close all;
clc;
image_name='garden1.jpg';
% Show Original Image
figure;
hold on;
imshow(image_name);
title({['Original Image']});
axis equal;
print(gcf,'-djpeg' ,strcat('original_image.jpeg'),'-r400')
% Choose Standard Deviation (sigma value = 1)
sigma = 1;
blur_image=Gau_blur(sigma,image_name);
figure;
title({['Image after Gaussian Blur'];
    ['¦Ò=',num2str(sigma)]});
hold on;
imshow(blur_image);
axis equal;
print(gcf,'-djpeg' ,strcat('Question4_sigma_',num2str(sigma),'.jpeg'),'-r400')

% Choose Standard Deviation (sigma value = 3)
sigma = 3;
blur_image=Gau_blur(sigma,image_name);
figure;
title({['Image after Gaussian Blur'];
    ['¦Ò=',num2str(sigma)]});
hold on;
imshow(blur_image);
axis equal;
print(gcf,'-djpeg' ,strcat('Question4_sigma_',num2str(sigma),'.jpeg'),'-r400')

% Choose Standard Deviation (sigma value = 5)
sigma = 5;
blur_image=Gau_blur(sigma,image_name);
figure;
title({['Image after Gaussian Blur'];
    ['¦Ò=',num2str(sigma)]});
hold on;
imshow(blur_image);
axis equal;
print(gcf,'-djpeg' ,strcat('Question4_sigma_',num2str(sigma),'.jpeg'),'-r400')
end

function blur_image=Gau_blur(sigma,image_name)
% Read in the Image
Image = imread(image_name);
% Change Format to Double
Img = double(Image);
% Calculate the half value of Kernel size
half_size=floor(3*sigma);
% Form the Gaussian Kernel
[x,y]=meshgrid(-half_size:half_size,-half_size:half_size);
Exp_index = -(x.^2+y.^2)/(2*sigma*sigma);
Kernel= exp(Exp_index)/(2*pi*sigma*sigma);
% Initialize
Output=zeros(size(Img));
% Extend the Image Border Pixels with Zeros
Img = padarray(Img,[half_size half_size]);
% Do Convolution for Each Color
for color=1:size(Img,3)
    for i = 1:size(Output,1)
        for j =1:size(Output,2)
            temp = Img(i:i+2*half_size,j:j+2*half_size,color).*Kernel;
            Output(i,j,color)=sum(temp(:));
        end
    end
end
% Image without Noise after Gaussian blur
blur_image = uint8(Output);
end