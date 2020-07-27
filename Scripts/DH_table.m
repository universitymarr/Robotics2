clear all
%% Define symbolic variables
syms alpha a d C D E F theta
syms q1 q2 q3 q4 q5 q6 real
syms  d1 d2 d3 a2 a3 d4 L1 L2 L3 L4 r2 real

%% number of joints
N=6;

%% Insert DH table of parameters of SCARA

%% ALPHA A D THETA
DHTABLE = [-pi/2 0 L1 q1;
           0 L2 0 q2;
           0 L3 0 q3;
           -pi/2 0 0 q4;
           0 L4 0 q5;
           0 0 0 q6]
    
%% Build the general Denavit-Hartenberg trasformation matrix
TDH = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
       sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
          0             sin(alpha)             cos(alpha)            d;
          0               0                      0                   1];
      
%% Build transformation matrices for each link

% empty cell array, we will put homogenous matrices
A = cell(1,N);  

% For every row in 'DHTABLE' we substitute the right value inside
% the general DH matrix
for i = 1:N
    alpha = DHTABLE(i,1);
    a = DHTABLE(i,2);
    d = DHTABLE(i,3);
    theta = DHTABLE(i,4);
    A{i} = subs(TDH);
end

% we compute the final result matrix
% homogenous matrix of last frame respect to the first frame
tot=A{1};
for i = 2:N
    ris=A{i};
    tot=tot*ris;
end

%% Compute the point of last frame respect to the first frame
p=simplify(tot(1:3,4))
pause;

%% Compute matrices R from homogenous A
R = cell(1,N);  

for i = 1:N
    R{i}=A{i}(1:3,1:3);
    fprintf('R{%d} \n',i);
    R{i}
    pause;
end
