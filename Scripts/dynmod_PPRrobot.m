%% dynamic model of PPR planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 m3 real
syms d1 d2 d3 real
syms L1 L2 L3 real
syms I1xx I1yy I1zz real
syms I2xx I2yy I2zz real  
syms I3xx I3yy I3zz real  
syms q1 q2 q3 real
syms dq1 dq2 dq3 real 
syms ddq1 ddq2 ddq3 real 
syms u1 u2  u3 g0 real

disp('**** dynamic model of PPR planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

disp(' ')
disp('*kinetic energy of link 1*')

%% compute kinetic energy of joint 1
pc1 = [0 q1-d1 0]';
vc1 = diff(pc1,q1)*dq1;

T1 = simplify((1/2)*m1*vc1'*vc1)
pause

disp('*kinetic energy of link 2*')

%% compute kinetic energy of joint 2
pc2 = [q2-d2 q1 0]';
vc2 = diff(pc2,q1)*dq1 + diff(pc2,q2)*dq2;

T2 = simplify((1/2)*m2*vc2'*vc2)
pause

disp('*kinetic energy of link 3*')

%% compute linear part of kinetic energy of joint 3
pc3 = [q2+d3*cos(q3) q1+d3*sin(q3) 0]';
vc3 = simplify(diff(pc3, q1)*dq1+diff(pc3,q2)*dq2 + diff(pc3,q3)*dq3);
Tl3 = simplify((1/2)*m3*vc3'*vc3);

%% compute angular part of kinetic energy of joint 3
om3 = [0 0 dq3]';
Ta3 = simplify((1/2)*om3'*diag([I3xx I3yy I3zz])*om3);

T3 = simplify(Tl3+Ta3)
pause

T=simplify(T1+T2+T3);

T = collect(T,dq1^2);
T = collect(T,dq2^2);
T = collect(T,dq3^2)
pause

%% compute inertia matrix 

M(1,1)=diff(T,dq1,2);
M(2,2)=diff(T,dq2,2);
M(3,3)=diff(T,dq3,2);

TempB12=diff(T,dq1);
M(1,2)=diff(TempB12,dq2);
M(2,1)=M(1,2);
M(1,3)=diff(TempB12,dq3);
M(3,1)=M(1,3);

TempB23 = diff(T,dq2);
M(2,3)=diff(TempB23,dq3);
M(3,2)=M(2,3);

M = simplify(M)
pause

%% parametrization
syms a1 a2 a3 a4 real

M=[a1 0 a4*cos(q3);
    0 a2 -a4*sin(q3);
    a4*cos(q3) -a4*sin(q3) a3;]

%% compute the Christoffel Matrix
q=[q1;q2;q3];

M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1))
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2))
M3=M(:,3);
C3=(1/2)*(jacobian(M3,q)+jacobian(M3,q)'-diff(M,q3))

pause

disp('*robot centrifugal terms*')

dq=[dq1;dq2;dq3];
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c3=dq'*C3*dq;
c=[c1;c2;c3]

pause 

disp('*potential energy of link 1*')
%% gravity vector
g= [0;g0;0];

%% compute the potential energy of link 1
U1=0

disp('*potential energy of link 2*')

%% compute the potential energy of link 2
U2=0;

disp('* total potential energy*')

%% total potential energy
U=U1+U2

%% G
G= jacobian(U,q)'

%% show complete dynamic equations
M*[ddq1; ddq2; ddq3]+c+G==[u1 u2 u3]'

syms V0 alpha real

qdd =[((L3*sin(q3)*dq3^2 + V0*cos(alpha))*(L3^2*sin(q3)^2 + 1))/(L3^2 + 1) - (L3^2*cos(q3)*sin(q3)*(- L3*cos(q3)*dq3^2 + V0*sin(alpha)))/(L3^2 + 1);
 (L3^2*cos(q3)*sin(q3)*(L3*sin(q3)*dq3^2 + V0*cos(alpha)))/(L3^2 + 1) - ((- L3*cos(q3)*dq3^2 + V0*sin(alpha))*(L3^2 - L3^2*sin(q3)^2 + 1))/(L3^2 + 1);
 L3*V0*cos(alpha-q3)]
  

disp('***end***')
