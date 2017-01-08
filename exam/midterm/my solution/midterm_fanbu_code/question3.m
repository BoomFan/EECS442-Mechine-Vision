function main()
clear all;
close all;
clc;
[points1 points2] = readPoints('4points.txt')
H=DLT_H(points1,points2)
end
function H=DLT_H(x1,x2)
[n1, ~]=size(x1);
[n2, ~]=size(x2);
if n1~=n2
    error=char('x1 and x2 does not match!')
   return
else
    n=n1;
end

%Build the matrix A such that,
% A_i=[-x,-y,-1,0,0,0,ux,uy,u;
%     0,0,0,-x,-y,-1,vx,vy,v]
A=zeros(2*n,9);
for i = 1:n
    A(2*i-1,1:3)=[-x1(i,1),-x1(i,2),-1];
    A(2*i-1,7:9)=[x2(i,1)*x1(i,1),x2(i,1)*x1(i,2),x2(i,1)];
    A(2*i,4:6)=[-x1(i,1),-x1(i,2),-1];
    A(2*i,7:9)=[x2(i,2)*x1(i,1),x2(i,2)*x1(i,2),x2(i,2)];
end
h=null(A);
H=zeros(3,3);
H(1,:)=h(1:3);
H(2,:)=h(4:6);
H(3,:)=h(7:9);
H=H/H(3,3);

% % If more than 4 points, we can use SVD to solve H
% [u,s,v] = svd(A,0);
% vv=v(:,9);
% for i=1:3
% H(1,i)=vv(i);
% end
% for i=1:3
% H(2,i)=vv(i+3);
% end
% for i=1:3
% H(3,i)=vv(i+6);
% end
% 
% % let rank(F)=2
% [u,s,v] = svd(H);
% H = H - u(:,3)*s(3,3)*transpose(v(:,3));
% H = H/H(3,3);
end
