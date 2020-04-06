%% dynamic model of RP planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 real
syms dc1 dc2 dc3 real
syms l1 l2 l3 real
syms I1xx I1yy I1zz real  %symbolic variables explicitly defined as real
syms I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms q1 q2 real
syms dq1 dq2 real 
syms ddq1 ddq2 real 
syms u1 u2 g0 real

disp('**** dynamic model of RP planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

disp(' ')
disp('*kinetic energy of link 1*')
%% compute linear part of linear kinetic energy of joint 1
pc1 = [dc1*cos(q1); dc1*sin(q1);0];
vc1=diff(pc1,q1)*dq1;
Tl1 = 1/2*m1*vc1'*vc1;



%% compute angular part of linear kinetic energy of joint 1
om1 = [0 0 dq1]';
Ta1 = 1/2*om1'*diag([I1xx I1yy I1zz])*om1;

T1 = Tl1+Ta1
pause

disp('*kinetic energy of link 2')

%% compute the linear part of kinetic energy of joint 2
pc2=[l1*cos(q1)-(q2-dc2)*sin(q1); l1*sin(q1)+(q2-dc2)*cos(q1); 0];
vc2=diff(pc2,q1)*dq1+diff(pc2,q2)*dq2;
Tl2= (1/2)*m2*vc2'*vc2;

%% compute the angular part of kinetic energy of joint 2
om2 = [0 0 dq1]';
Ta2 = 1/2*om2'*diag([I2xx I2yy I2zz])*om2;

T2=simplify(Tl2+Ta2)
pause

disp('***robot kinetic energy***')

%% total kinetic energy
T=T1+T2;

disp('*simplifying*')

T=simplify(T1+T2)
pause

% collect in base of term that you pass it in this case, collect terms that has dq1^2 and 
% do the same for dq2^2 because you need them to compute M
T=collect(T,dq1^2);
T=collect(T,dq2^2);

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
