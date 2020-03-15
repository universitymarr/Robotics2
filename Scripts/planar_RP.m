%% RP planar
clear all

syms q1 q2
syms qd2 qd3

vec =[q1 q2];

%% define point 
p=[q2*cos(q1);
    q2*sin(q1)];

%% Analytic Jacobian
J=simplify(jacobian(p,vec))
pause

%% define velocity in Cartesian Space
pd= []; 

%% define acceleration in Cartesian Space
pdd=[];

%% time derivative of the Jacobian PPR
Jd_RP =[cos(q1)*(qd2)-sin(q1)*qd1*q2;
        sin(q1)*(qd2)+cos(q1)*qd1*q2];
