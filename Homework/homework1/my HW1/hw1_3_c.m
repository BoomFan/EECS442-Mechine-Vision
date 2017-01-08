clear all;
clc;
syms d x y z a b c

H=[1.520,-1.902,1;
    3.3,23.49,3;
    1,3,1]
x1=[1;3;1]
x2=[3;1;1]
xx1=H*x1
xx2=H*x2

X=[x;-x+4;1]
XX=H*X
pretty(XX)


[a,b]=solve('a*((1711*0)/500 - 826/125)+b=2424/25 - (2019*0)/100','a*((1711*1)/500 - 826/125)+b=2424/25 - (2019*1)/100')

format short
a
b

detH=det(H)
HH=transpose(H^-1)
new_H=detH*HH