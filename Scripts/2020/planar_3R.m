%% 3R planar
clear all
clc
syms q1 q2 q3 real
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms l1 l2 l3 real
syms l real

vec =[q1 q2 q3];
%% define length of joints
l1 = 1;
l2= 1;
l3 =1;

%% define point 
p=[l1*cos(q1)+l2*cos(q1+q2)+l3*cos(q1+q2+q3);
    l1*sin(q1)+l2*sin(q1+q2)+l3*sin(q1+q2+q3);];
    
%% Analytic Jacobian
J=simplify(jacobian(p,vec))
pause

%% define velocity in Cartesian Space
pd= [1; 1]; 

%% define acceleration in Cartesian Space
pdd=[4;2];

%% time derivative of the Jacobian 3R
dJ = simplify(diff(J,q1)*dq1 + diff(J,q2)*dq2 + diff(J,q3)*dq3)
pause

ddJ = simplify(diff(dJ,q1)*dq1 + diff(dJ,q2)*dq2 + diff(dJ,q3)*dq3 + diff(dJ,dq1)*ddq1 + diff(dJ,dq2)*ddq2 + diff(dJ,dq3)*ddq3)
