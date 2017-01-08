function main
clear all;
close all;
clc;
% load data
set_number=1;
[x1, x2] = readTextFiles(strcat('set',num2str(set_number)));   % default setting is to open set1 data
image1 = imread(strcat('set',num2str(set_number),'/image1.jpg'));
image2 = imread(strcat('set',num2str(set_number),'/image2.jpg'));

% Linear least squares version
F_lin = cal_F(x1,x2);
% Normalized version
T1=cal_T(x1);
T2=cal_T(x2);
xt1=T1*x1;
xt2=T2*x2;
hold on
plot(x1(1,:),x1(2,:),'bo')
figure
plot(xt1(1,:),xt1(2,:),'ro')
F_temp=cal_F(xt1,xt2);
F_normal=transpose(T1)*F_temp*(T2);
% Epipolar line and error
[L_linear_1,L_linear_2,error_lin_1,error_lin_2] = epi_line_error(x1,x2,F_lin);
[L_norm_1,L_norm_2,error_norm_1,error_norm_2] = epi_line_error(x1,x2,F_normal);
%visualization
figure;
hold on;
draw_pic_linear(x1,x2,L_linear_1,L_linear_2,error_lin_1,error_lin_2,image1,image2,set_number)
figure;
hold on;
draw_pic_normal(x1,x2,L_norm_1,L_norm_2,error_norm_1,error_norm_2,image1,image2,set_number)
end

function draw_pic_normal(x1,x2,L1,L2,error1,error2,image1,image2,set_number)
[~, n]=size(x1);
line_len=15;
h_title=suptitle({['Fundamental Matrix'],
    ['Result for Set',num2str(set_number)]});
subplot(1,2,1)
hold on;
h_title=title({['Image1 normalized version'];
    ['Average distance=',num2str(error1)]});
imshow(image1);
plot(x1(1,:),x1(2,:),'ro');
for i = 1:n
    if L1(2,i)==0
        p1 = [-L1(3,i)/L1(1,i),x1(2,i)-line_len];
        p2 = [-L1(3,i)/L1(1,i),x1(2,i)+line_len];
    else
        p1 = [x1(1,i)-line_len,x1(1,i)+line_len];
        p2 = [-(L1(1,i)*p1(1,1)+L1(3,i))/L1(2,i), -(L1(1,i)*p1(1,2)+L1(3,i))/L1(2,i)]; 
    end
        plot(p1,p2,'b');
end
% Plot image2
subplot(1,2,2)
hold on;
h_title=title({['Image2 normalized version'];
    ['Average distance=',num2str(error2)]});
imshow(image2);
plot(x2(1,:),x2(2,:),'ro');
for i = 1:n
    if L2(2,i)==0
        p1 = [-L2(3,i)/L2(1,i),x2(2,i)-line_len];
        p2 = [-L2(3,i)/L2(1,i),x2(2,i)+line_len];
    else
        p1 = [x2(1,i)-line_len,x2(1,i)+line_len];
        p2 = [-(L2(1,i)*p1(1,1)+L2(3,i))/L2(2,i), -(L2(1,i)*p1(1,2)+L2(3,i))/L2(2,i)]; 
    end
        plot(p1,p2,'b');
end
print(gcf,'-djpeg' ,strcat('HW3_2_1_normalized_set',num2str(set_number),'.jpeg'),'-r400')

end

function draw_pic_linear(x1,x2,L1,L2,error1,error2,image1,image2,set_number)
[~, n]=size(x1);
line_len=15;
% Plot image1
h_title=suptitle({['Fundamental Matrix'],
    ['Result for Set',num2str(set_number)]});
subplot(1,2,1)
hold on;
h_title=title({['Image1 linear least square version'];
    ['Average distance=',num2str(error1)]});
imshow(image1);
plot(x1(1,:),x1(2,:),'ro');
for i = 1:n
    if L1(2,i)==0
        p1 = [-L1(3,i)/L1(1,i),x1(2,i)-line_len];
        p2 = [-L1(3,i)/L1(1,i),x1(2,i)+line_len];
    else
        p1 = [x1(1,i)-line_len,x1(1,i)+line_len];
        p2 = [-(L1(1,i)*p1(1,1)+L1(3,i))/L1(2,i), -(L1(1,i)*p1(1,2)+L1(3,i))/L1(2,i)]; 
    end
        plot(p1,p2,'b');
end
% Plot image2
subplot(1,2,2)
hold on;
h_title=title({['Image2 linear least square version'];
    ['Average distance=',num2str(error2)]});
imshow(image2);
plot(x2(1,:),x2(2,:),'ro');
for i = 1:n
    if L2(2,i)==0
        p1 = [-L2(3,i)/L2(1,i),x2(2,i)-line_len];
        p2 = [-L2(3,i)/L2(1,i),x2(2,i)+line_len];
    else
        p1 = [x2(1,i)-line_len,x2(1,i)+line_len];
        p2 = [-(L2(1,i)*p1(1,1)+L2(3,i))/L2(2,i), -(L2(1,i)*p1(1,2)+L2(3,i))/L2(2,i)]; 
    end
        plot(p1,p2,'b');
end
print(gcf,'-djpeg' ,strcat('HW3_2_1_LinearLS_set',num2str(set_number),'.jpeg'),'-r400')
end


function [L1,L2,error_1,error_2]=epi_line_error(x1,x2,F)
[~, n]=size(x1);
L1 = F*x2;
L2 = transpose(F)*x1;
% distance=|ax+by+c|/sqrt(a^2+b^2)
err1=sum(L1.*x1);   % calculate ax+by+c
den1=sqrt((L1(1,:).^2)+L1(2,:).^2); % calculate denominator
dist1=err1./den1;   % calculate each distance
err2=sum(L2.*x2);
den2=sqrt((L2(1,:).^2)+L2(2,:).^2);
dist2=err2./den2;
error_1=sum(abs(dist1))/n;
error_2=sum(abs(dist2))/n;
end

%% Calculate Transformation Matrix
function T=cal_T(x)
[~, n]=size(x);
x_bar=sum(x(1,:))/n;
y_bar=sum(x(2,:))/n;
i=1;
num=sqrt((x(1,i)-x_bar)^2+(x(2,i)-y_bar)^2);
den=n*sqrt(2);
d=num/den;
if n>=2
    for i=2:n
        num=sqrt((x(1,i)-x_bar)^2+(x(2,i)-y_bar)^2);
        den=n*sqrt(2);
        d=d+num/den;
    end
else
end
T=[1/d,0,-x_bar/d;
    0,1/d,-y_bar/d;
    0,0,1];

end

%% Calculate Fundamental Matrix
function F=cal_F(x1,x2)
[~, n1]=size(x1);
[~, n2]=size(x2);
if n1~=n2
    error=char('x1 and x2 does not match!')
   return
else
    n=n1;
end

%Build the matrix A
for i = 1:n
    xx1 = x1(:,i);
    xx2 = x2(:,i);
    xx=xx2*transpose(xx1);
    for j=1:9
        A(i,j)=xx(j);
    end
end

%SVD
[u,s,v] = svd(A,0);
vv=v(:,9);
for i=1:3
F(1,i)=vv(i);
end
for i=1:3
F(2,i)=vv(i+3);
end
for i=1:3
F(3,i)=vv(i+6);
end

% let rank(F)=2
[u,s,v] = svd(F);
F = F - u(:,3)*s(3,3)*transpose(v(:,3));
end