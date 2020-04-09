%% dynamic model of RPP planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 m3 real
syms d1 d2 d3 real
syms r2 a2 real
syms I1xx I1yy I1zz real
syms I2xx I2yy I2zz real 
syms I3xx I3yy I3zz real
syms q1 q2 q3 real 
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms u1 u2 u3 g0 real

disp('**** dynamic model of RPP planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

disp(' ')
disp('*kinetic energy of link 1*')

%% compute  kinetic energy of joint 1
om1 = [0 0 dq1]';
T1=(1/2)*om1'*diag([I1xx I1yy I1zz])*om1

pause

disp('*kinetic energy of link 2*')

%% compute linear part kinetic energy of joint 2
pc2 = [r2*q1 q2 0]';
vc2 = diff(pc2,q1)*dq1+diff(pc2,q2)*dq2;
Tl2 = (1/2)*m2*vc2'*vc2;

%% compute  kinetic energy of joint 2
om2 = [0 0 dq1]';
Ta2=(1/2)*om2'*diag([I2xx I2yy I2zz])*om2;

T2 = simplify(Tl2+Ta2)
pause

disp('*kinetic energy of link 3*')

%% compute linear part kinetic energy of joint 3
tot_v = (a2^2+q3^2)*dq1^2+dq2^2+dq3^2+2*a2*dq1*dq3;
Tl3 = (1/2)*m3*tot_v;

%% compute  kinetic energy of joint 2
om3 = [0 0 dq1]';
Ta3=(1/2)*om3'*diag([I3xx I3yy I3zz])*om3;

T3 = simplify(Tl3+Ta3)
pause

T=simplify(T1+T2+T3);

T=collect(T,dq1^2);
T=collect(T,dq2^2);
T=collect(T,dq3^2)
pause

disp('***robot inertia matrix***')

%% compute robot matrix M, inertia matrix

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

M=simplify(M)

%% parametrization
syms b1 b2 b3 b4 real

Mp=[b1+b3*q3^2 0 b4;
   0 b2 0;
   b4 0 b3;]

%% compute the Christoffel Matrix
q=[q1;q2;q3];

M1p=Mp(:,1);
C1p=(1/2)*(jacobian(M1p,q)+jacobian(M1p,q)'-diff(Mp,q1))
M2p=Mp(:,2);
C2p=(1/2)*(jacobian(M2p,q)+jacobian(M2p,q)'-diff(Mp,q2))
M3p=Mp(:,3);
C3p=(1/2)*(jacobian(M3p,q)+jacobian(M3p,q)'-diff(Mp,q3))

pause

disp('***robot centrifugal terms (no Coriolis here!)***')

dq=[dq1;dq2;dq3];
c1p=dq'*C1p*dq;
c2p=dq'*C2p*dq;
c3p=dq'*C3p*dq;
cp=[c1p;c2p;c3p]

pause 

disp('*potential energy of link 1*')

%% compute the potential energy of link 1
U1=0;

disp('*potential energy of link 2*')

%% vector gravity acceleration
g=[0;-g0;0];

%% compute the potential energy of link 2
U2=-m2*g*q2;

disp('*potential energy of link 3*')

%% compute the potential energy of link 2
U3=-m3*g*q2;

U=U1+U2+U3
pause

G= jacobian(U,q)
pause