%% 2R planar
clear all

syms q1 q2 real
syms dq1 dq2 real
syms ddq1 ddq2 real
syms l1 l2 l3 real


vec =[q1 q2];
%% define length of joints
l1 = 1;
l2= 1;
%% define point 
p=[l1*cos(q1)+l2*cos(q1+q2);
    l1*sin(q1)+l2*sin(q1+q2);];

%% Analytic Jacobian
J=simplify(jacobian(p,vec))
pause

%% define velocity in Cartesian Space
pd= [1; 1]; 

%% define acceleration in Cartesian Space
pdd=[4;2];

%% time derivative of the Jacobian 3R
dJ=simplify(diff(J,q1)*dq1 + diff(J,q2)*dq2)

pause

ddJ = simplify(diff(dJ,q1)*dq1 + diff(dJ,q2)*dq2 + diff(dJ,dq1)*ddq1 + diff(dJ,dq2)*ddq2)


