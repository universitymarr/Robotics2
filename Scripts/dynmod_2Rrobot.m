%% dynamic model of 2R planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 real
syms L1 L2 real
syms d1 d2 real
syms I1xx I1yy I1zz real  %symbolic variables explicitly defined as real
syms I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms q1 q2 real
syms dq1 dq2 real
syms ddq1 ddq2 real
syms u1 u2 g0 real

disp('**** dynamic model of 2R planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

disp(' ')
disp('*kinetic energy of link 1 - linear part*')

%% compute linear part of kinetic energy of joint 1
pc1 =[d1*cos(q1) d1*sin(q1) 0]';
vc1 = diff(pc1,q1)*dq1;
Tl1 = (1/2)*m1*vc1'*vc1;

disp('*kinetic energy of link 1 - angular part*')

%% compute the angular part of kinetic energy of joint 2
om1 = [0 0 dq1]';
Ta1 = (1/2)*om1'*diag([I1xx I1yy I1zz])*om1;

T1= simplify(Tl1+Ta1)
pause

disp('*kinetic energy of link 2 - linear part*')

%% compute the linear part of kinetic energy of joint 2
pc2 = [L1*cos(q1)+d2*cos(q2) L1*sin(q1)+d2*sin(q2) 0]';
vc2 = diff(pc2,q1)*dq1+diff(pc2,q2)*dq2;
T2l = (1/2)*m2*vc2'*vc2;

disp('*kinetic energy of link 2 - angular part*')

%% compute the angular part of kinetic energy of joint 2
om2=[0 0 dq2]';
T2a=(1/2)*om2'*diag([I2xx I2yy I2zz])*om2;

%% kinetic energy of joint 2
T2=simplify(T2l+T2a)
pause

disp('***robot kinetic energy***')
    
%% total kinetic energy
T=T1+T2;

disp('*simplifying*')
    
T=simplify(T1+T2);

% collect in base of term that you pass it in this case, collect terms that has dq1^2 and 
% do the same for dq2^2 because you need them to compute M
T=collect(T,dq1^2);
T=collect(T,dq2^2)
pause

disp('***robot inertia matrix***')

%% compute robot matrix M, inertia matrix

% we extract the element M(1,1) from the total kinetic energy by doing twice the derivative
% respect qd1 because this will be the factor that multiplyin the square of velocity  of joint 1
% we do the same also for joint 2
% the factor 1/2 disappear automatically in this differentation
M(1,1)=diff(T,dq1,2);
M(2,2)=diff(T,dq2,2);
TempB12=diff(T,dq1);
M(1,2)=diff(TempB12,dq2);
M(2,1)=M(1,2);
M=simplify(M)

pause

%% parametrization 
syms a1 a2 a3 a4 real
M=[ a1 a3*cos(q1-q2);
    a3*cos(q1-q2) a2;]
disp('*Christoffel matrices*')

%% compute the Christoffel Matrix
q=[q1;q2];

% consider the first column of matrix M, so you have M1
% do the same for M2
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1))
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2))

pause

disp('***robot centrifugal terms (no Coriolis here!)***')

dq=[dq1;dq2];
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c=[c1;c2]

pause 

disp('*potential energy of link 1*')

%% compute the potential energy of link 1
U1=0

pause

disp('*potential energy of link 2*')

%% vector gravity acceleration
g=[0;-g0;0];

%% compute the potential energy of link 2
U2=0
pause

disp('***robot potential energy (due to gravity)***')

%%  total potential energy
U=U1+U2

pause

disp('***robot gravity term***')

%% compute G
G=jacobian(U,q)'

pause

disp('***complete dynamic equations***')

%% show complete dynamic equations
M*[ddq1; ddq2]+c+G==[u1 u2]'

disp('***end***')

% end
