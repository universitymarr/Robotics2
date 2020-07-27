%% 4R planar
clear all
clc
syms q1 q2 q3 q4 real
syms dq1 dq2 dq3 dq4 real
syms l1 l2 l3 l4 real


vec =[q1 q2 q3 q4];
%% define length of joints
% l1 = 0.2;
% l2 = 0.2;
% l3 = 0.2;
% l4 = 0.2;

%% define point 
p=[l1*cos(q1)+l2*cos(q1+q2)+l3*cos(q1+q2+q3)+ l4*cos(q1+q2+q3+q4);
    l1*sin(q1)+l2*sin(q1+q2)+l3*sin(q1+q2+q3)+l4*sin(q1+q2+q3+q4);];

%% Analytic Jacobian
J=simplify(jacobian(p,vec))

pause

%% define velocity in Cartesian Space
pd= [];

%% define acceleration in Cartesian Space
pdd=[5;2];

%% time derivative of the Jacobian 4R
dJ = diff(J,q1)*dq1 + diff(J,q2)*dq2 + diff(J,q3)*dq3 + diff(J,q4)*dq4
