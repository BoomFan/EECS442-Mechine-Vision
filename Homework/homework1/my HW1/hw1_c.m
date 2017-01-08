clear all;
clc;
syms d x y z a b c
T=[1,0,0,-d/sqrt(2);
    0,1,0,0;
    0,0,1,-d/sqrt(2);
    0,0,0,1]
theta=3*pi/4
%���������y�����ת�������������Ϊ��ʵ����y�������򣬶���x�ᵽz������ֶ�����
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
