%% dynamic model of Boulton-Watt under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m real
syms d1 d2 d3 real
syms L h real
syms Isxx Isyy Iszz real  %symbolic variables explicitly defined as real
syms I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms theta phi real
syms dtheta dphi real 
syms ddq1 ddq2 real 
syms u1 u2 g0 real

disp(' ')
disp('*kinetic energy of link 1*')
%% compute linear part of linear kinetic energy of joint 1
pct = [0;0;h];
vct=diff(pct,theta)*dtheta;
Tlt = simplify(1/2*m*vct'*vct);

%% compute angular part of linear kinetic energy of joint 1
omt = [0 0 dtheta]';
Tat = simplify(1/2*omt'*diag([Isxx Isyy Iszz])*omt);

Tt = simplify(Tlt+Tat)
pause

disp('*kinetic energy of link 2')

%% compute the linear part of kinetic energy of joint 2
%pcp=[L*cos(phi)*-theta L*sin(phi) 0]';
%vcp=diff(pcp,theta)*dtheta+diff(pcp,phi)*dphi;
%vcp = dtheta*sin(phi)+phi;
Tlp= simplify((1/2)*m*L^2*(dphi^2+dtheta^2*sin(phi)^2));

%% compute the angular part of kinetic energy of joint 2
omp = [0 0 0]';
Tap = simplify(1/2*omp'*diag([Isxx Isyy Iszz])*omp);

Tp= simplify(Tlp+Tap)
pause

disp('***robot kinetic energy***')

%% total kinetic energy
T=simplify(Tt+2*Tp);
T=collect(T,dtheta^2);
T=collect(T,dphi^2);

%% inertia matrix

M(1,1)=diff(T,dtheta,2);
M(2,2)=diff(T,dphi,2);

TempB12=diff(T,dtheta);
M(1,2)=diff(TempB12,dphi);
M(2,1)=M(1,2);
M=simplify(M)
pause
%% compute the Christoffel Matrix
q=[theta;phi];

% consider the first column of matrix M, so you have M1
% do the same for M2
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,theta))
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,phi))

pause

disp('***robot centrifugal terms (no Coriolis here!)***')

dq=[dtheta;dphi];
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c=[c1;c2]

pause 

%% gavity term
g=-g0;


U = 2*m*g*L*cos(phi);

G=jacobian(U,q)'
pause


