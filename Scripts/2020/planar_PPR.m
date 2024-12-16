%% PPR planar
clear all

syms q1 q2 q3
syms qd3

vec =[q1 q2 q3];
%% define length of joints
syms l

%% define point 
p=[q1+l*cos(q3);
    q2+l*sin(q3)];

%% Analytic Jacobian
J=simplify(jacobian(p,vec))
pause

%% define velocity in Cartesian Space
pd= [-1;1]; 

%% define acceleration in Cartesian Space
pdd=[];

%% time derivative of the Jacobian PPR
Jd_PPR =  [0 0 -l*cos(q3)*qd3;
           0 0 -l*sin(q3)*qd3];
