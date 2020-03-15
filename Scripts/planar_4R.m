%% 4R planar
clear all

syms q1 q2 q3 q4
syms qd1 qd2 qd3 qd4
syms l1 l2 l3 l4


vec =[q1 q2 q3 q4];
%% define length of joints
% l1 = 1;
% l2= 1;
% l3 =1;
% l4 = 1

%% define point 
p=[l1*cos(q1)+l2*cos(q1+q2)+l3*cos(q1+q2+q3)+ l4*cos(q1+q2+q3+q4);
    l1*sin(q1)+l2*sin(q1+q2)+l3*sin(q1+q2+q3)+l4*sin(q1+q2+q3+q4);];

%% Analytic Jacobian
J=simplify(jacobian(p,vec))
pause

%% define velocity in Cartesian Space
pd= [];

%% define acceleration in Cartesian Space
pdd=[4;2];

%% time derivative of the Jacobian 4R
Jd_4R=[-l1*cos(q1)*qd1-l2*cos(q1+q2)*(qd1+qd2)-l3*cos(q1+q2+q3)*(qd1+qd1+qd3)-l4*cos(q1+q2+q3+q4)*(qd1+qd2+qd3+qd4) -l2*cos(q1+q2)*(qd1+qd2)-l3*cos(q1+q2+q3)*(qd1+qd2+qd3)-l4*cos(q1+q2+q3+qd4)*(qd1+qd2+qd3+qd4) -l3*cos(q1+q2+q3)*(qd1+qd2+qd3)-l4*cos(q1+q2+q3+q4)*(qd1+qd2+qd3+qd4) -l4*cos(q1+q2+q3+q4)*(qd1+qd2+qd3+qd4);
      -l1*sin(q1)*qd1-l2*sin(q1+q2)*(qd1+qd2)-l3*sin(q1+q2+q3)*(qd1+qd2+qd3)-l4*sin(q1+q2+q3+q4)*(qd1+qd2+qd3+qd4) -l2*sin(q1+q2)*(qd1+qd2)-l3*sin(q1+q2+q3)*(qd1+qd2+qd3)-l4*sin(q1+q2+q3+q4)*(qd1+qd2+qd3+qd4) -l3*sin(q1+q2+q3)*(qd1+qd2+qd3)-l4*sin(q1+q2+q3+q4)*(qd1+qd2+qd3+qd4) l4*sin(q1+q2+q3+q4)*(qd1+qd2+qd3+qd4)];

