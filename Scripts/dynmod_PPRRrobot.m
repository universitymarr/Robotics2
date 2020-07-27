%% dynamic model of PPRR planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 m3 m4 real
syms d1 d2 d3 d4 real
syms L1 L2 L3 L4 real
syms I1xx I1yy I1zz real  %symbolic variables explicitly defined as real
syms I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms I3xx I3yy I3zz real  %symbolic variables explicitly defined as real
syms I4xx I4yy I4zz real  %symbolic variables explicitly defined as real
syms q1 q2 q3 q4 real
syms dq1 dq2 dq3 dq4 real 
syms ddq1 ddq2 ddq3 ddq4 real 
syms u1 u2 u3 u4 g0 real

disp('**** dynamic model of RRRR planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

disp(' ')
disp('*kinetic energy of link 1*')

%% compute kinetic energy of joint 1
pc1 = [-(d1-q1); 0;0];
vc1=diff(pc1,q1)*dq1;
T1 = 1/2*m1*vc1'*vc1
pause
disp('*kinetic energy of link 2')

%% compute  kinetic energy of joint 2
pc2 = [q1; q2-d2;0];
vc2=diff(pc2,q1)*dq1+diff(pc2,q2)*dq2;
T2 = 1/2*m2*vc2'*vc2
pause

disp('*kinetic energy of link 3')

%% compute linear part of linear kinetic energy of joint 3
pc3 = [q1+d3*cos(q3); q2+d3*sin(q3);0];
vc3=diff(pc3,q1)*dq1+diff(pc3,q2)*dq2++diff(pc3,q3)*dq3;
Tl3 = 1/2*m3*vc3'*vc3;


%% compute the angular part of kinetic energy of joint 3
om3 = [0 0 dq3]';
Ta3 = 1/2*om3'*diag([I3xx I3yy I3zz])*om3;

T3= simplify(Tl3+Ta3)
pause
disp('*kinetic energy of link 4')

%% compute linear part of linear kinetic energy of joint 4
pc4 = [q1+L3*cos(q3)+d4*cos(q4+q3); q2+L3*sin(q3)+d4*sin(q4+q3);0];
vc4=diff(pc4,q1)*dq1+diff(pc4,q2)*dq2+diff(pc4,q3)*dq3+diff(pc4,q4)*dq4;
Tl4 = 1/2*m4*vc4'*vc4;

%% compute the angular part of kinetic energy of joint 3
om4 = [0 0 dq4+dq3]';
Ta4 = 1/2*om4'*diag([I4xx I4yy I4zz])*om4;

T4= simplify(Tl4+Ta4)
pause

T=simplify(T1+T2+T3+T4)
pause

T=collect(T,dq1^2);
T=collect(T,dq2^2);
T=collect(T,dq3^2);
T=collect(T,dq4^2);
disp('***robot inertia matrix***')

%% compute robot matrix M, inertia matrix

M(1,1)=diff(T,dq1,2);
M(2,2)=diff(T,dq2,2);
M(3,3)= diff(T,dq3,2);
M(4,4)=diff(T,dq4,2);

TempB1=diff(T,dq1);
M(1,2)=diff(TempB1,dq2);
M(2,1)=M(1,2);

M(1,3)=diff(TempB1,dq3);
M(3,1)=M(1,3);

M(1,4)=diff(TempB1,dq4);
M(4,1)=M(1,4);

TempB2=diff(T,dq2);
M(2,3) = diff(TempB2,dq3);
M(3,2)=M(2,3);

M(2,4)=diff(TempB2,dq4);
M(4,2)=M(2,4);

TempB3=diff(T,dq3);
M(3,4)=diff(TempB3,dq4);
M(4,3)=M(3,4);

M=simplify(M)
pause
disp('*parametrized*')

%% M in parametrized form

syms a1 a2 a3 a4 a5 a6 a7 real

M=[a1 0 -a6*sin(q3+q4)-a5*sin(q3) -a6*sin(q3+q4);
    0 a2 a6*cos(q3+q4)+a5*cos(q3) a6*cos(q3+q4);
    -a6*sin(q3+q4)-a5*sin(q3) a6*cos(q3+q4)+a5*cos(q3) a3+2*a6*L3*cos(q4) a4+a6*L3*cos(q4);
    -a6*sin(q3+q4) a6*cos(q3+q4) a4+a6*L3*cos(q4) a4;]

disp('*Christoffel matrices*')


%% compute the Christoffel Matrix
q=[q1;q2;q3;q4];


M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1))
pause

M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2))
pause

M3=M(:,3);
C3=(1/2)*(jacobian(M3,q)+jacobian(M3,q)'-diff(M,q3))
pause

M4=M(:,4);
C4=(1/2)*(jacobian(M4,q)+jacobian(M4,q)'-diff(M,q4))
pause

disp('***robot centrifugal terms (no Coriolis here!)***')

dq=[dq1;dq2;dq3;dq4];   
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c3=dq'*C3*dq;
c4=dq'*C4*dq;
c=[c1;c2;c3;c4]
pause

disp('*potential energy of link 1*')
%% vector gravity acceleration
g=[0;-g0;0];

%% compute the potential energy of link 1
U1=-m1*g'*pc1;

disp('*potential energy of link 2*')


%% compute the potential energy of link 2
U2=-m2*g'*pc2;

disp('*potential energy of link 3*')
%% compute the potential energy of link 3
U3=-m3*g'*pc3;

disp('*potential energy of link 4*')

%% compute the potential energy of link 4
U4=-m4*g'*pc4;

U=simplify(U1+U2+U3+U4)
pause

disp('***robot gravity term***')

G=jacobian(U,q)'
pause

%Gp = [0; a2*g0; g0*(a6*cos(q3)+a5*cos(q3+q4)); a5*g0]


disp('***complete dynamic equations***')

%% show complete dynamic equations
M*[ddq1; ddq2; ddq3; ddq4]+c+G==[u1 u2 u3 u4]'

ris = collect(M*[ddq1; ddq2; ddq3; ddq4]+c+G,a1);
ris = collect(ris,a2);
ris = collect(ris,a3);
ris = collect(ris,a4);
ris = collect(ris,a5);
ris = collect(ris,a6)

pause


