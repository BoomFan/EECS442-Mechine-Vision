clear all;
close all;
clc;
syms beta gamma
Ry=[cos(beta) 0 sin(beta);
    0 1 0;
    -sin(beta) 0 cos(beta)]
Rz=[cos(gamma) -sin(gamma) 0;
    sin(gamma) cos(gamma) 0;
    0 0 1;]
Rz*Ry

Ry*Rz
