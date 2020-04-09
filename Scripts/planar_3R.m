%% 3R planar
clear all

syms q1 q2 q3 real
syms dq1 dq2 dq3 real
%syms l1 l2 l3 real
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
dJ_3R=[-l1*cos(q1)*dq1-l2*cos(q1+q2)*(dq1+dq2)-l3*cos(q1+q2+q3)*(dq1+dq1+dq3) -l2*cos(q1+q2)*(dq1+dq2)-l3*cos(q1+q2+q3)*(dq1+dq2+dq3) -l3*cos(q1+q2+q3)*(dq1+dq2+dq3);
      -l1*sin(q1)*dq1-l2*sin(q1+q2)*(dq1+dq2)-l3*sin(q1+q2+q3)*(dq1+dq2+dq3) -l2*sin(q1+q2)*(dq1+dq2)-l3*sin(q1+q2+q3)*(dq1+dq2+dq3) -l3*sin(q1+q2+q3)*(dq1+dq2+dq3)];

