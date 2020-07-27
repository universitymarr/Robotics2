%% dynamic model of RPP planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 m3 real
syms d1 d2 d3 real
syms r2 a2 L1 real
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
pc2 = [r2*cos(q1) r2*sin(q1) L1+q2]';
vc2 = diff(pc2,q1)*dq1+diff(pc2,q2)*dq2;
Tl2 = (1/2)*m2*vc2'*vc2;

%% compute  kinetic energy of joint 2
om2 = [0 0 dq1]';
Ta2=(1/2)*om2'*diag([I2xx I2yy I2zz])*om2;

T2 = simplify(Tl2+Ta2)
pause

disp('*kinetic energy of link 3*')

%% compute linear part kinetic energy of joint 3
pc3 = [cos(q1)*(a2+q3) sin(q1)*(a2+q3) L1+q2+r2]'
%pc3 = [cos(q1)*(a2 + q3 + r2);  sin(q1)*(a2 + q3 + r2); -q2]


vc3 = diff(pc3,q1)*dq1+diff(pc3,q2)*dq2+diff(pc3,q3)*dq3;
Tl3 = (1/2)*m3*vc3'*vc3;

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
pause

%% parametrization
syms b1 b2 b3 b4 real

M=[b1+b3*(d3-q3)^2 0 b4;
    0 b2 0;
    b4 0 b3]

%% compute the Christoffel Matrix
q=[q1;q2;q3];

M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1))
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2))
M3=M(:,3);
C3=(1/2)*(jacobian(M3,q)+jacobian(M3,q)'-diff(M,q3))

pause

disp('***robot centrifugal terms (no Coriolis here!)***')

dq=[dq1;dq2;dq3];
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c3=dq'*C3*dq;
c=[c1;c2;c3]

pause 

disp('*potential energy of link 1*')

%% compute the potential energy of link 1
U1=0;

disp('*potential energy of link 2*')

%% vector gravity acceleration
g=[0;0;-g0];

%% compute the potential energy of link 2
U2=-m2*g'*pc2;

disp('*potential energy of link 3*')

%% compute the potential energy of link 2
U3=-m3*g'*pc3;

U=U1+U2+U3
pause

G= jacobian(U,q)'
pause