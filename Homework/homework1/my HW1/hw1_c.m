clear all;
clc;
syms d x y z a b c
T=[1,0,0,-d/sqrt(2);
    0,1,0,0;
    0,0,1,-d/sqrt(2);
    0,0,0,1]
theta=3*pi/4
%这里面相对y轴的旋转方向很暧昧，因为其实不是y轴正方向，而是x轴到z轴的右手定则方向
R=[cos(theta),0,sin(theta),0;
    0,1,0,0;
    -sin(theta),0,cos(theta),0;
    0,0,0,1]
P=R*T
pretty(P)

X=[x+a;y+b;z+c;1]
% X=[0;1;1;1]
% X=[0;y;z+1;1]
newX=P*X
pretty(newX)
