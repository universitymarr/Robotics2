%% dynamic model of 2R planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 m3 real
syms I1xx I1yy I1zz real
syms I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms I3xx I3yy I3zz real
syms d1 d2 d3 L1 L2 L3 real
syms q1 q2  q3 real
syms dq1 dq2 dq3  real 
syms ddq1 ddq2 ddq3 real 
syms u1 u2  u3 g0 real

disp('**** dynamic model of PRR planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause


%% link 1
pc1 = [0;0;q1-d1;];
vc1 = diff(pc1,q1)*dq1;
T1= (1/2)*m1*vc1'*vc1
pause

%% link 2
pc2 = [0;d2;q1;];
vc2 = diff(pc2,q1)*dq1 + diff(pc2,q2)*dq2;
Tl2= (1/2)*m2*vc2'*vc2;

om2= [0;dq2;0;];
Ta2= (1/2)*om2'*diag([I2xx I2yy I2zz])*om2;

T2=Ta2+Tl2
pause

%% link 3
pc3 = [d3*cos(q2)*cos(q3); L2+d3*sin(q3); q1-d3*cos(q3)*sin(q2);]
vc3 = diff(pc3,q1)*dq1 + diff(pc3,q2)*dq2 + diff(pc3,q3)*dq3;
Tl3 = (1/2)*m3*vc3'*vc3;

om3= [0;dq2;dq3;];
Ta3= (1/2)*om3'*diag([I3xx I3yy I3zz])*om3;

T3=Tl3+Ta3;

%% total kinetic energy
T=T1+T2+T3;

disp('*simplifying*')

T=simplify(T1+T2+T3);

% collect in base of term that you pass it in this case, collect terms that has dq1^2 and 
% do the same for dq2^2 because you need them to compute M
T=collect(T,dq1^2);
T=collect(T,dq2^2);
T=collect(T,dq3^2)
pause

disp('***robot inertia matrix***')

%% compute robot matrix M, inertia matrix

% we extract the element M(1,1) from the total kinetic energy by doing twice the derivative
% respect qd1 because this will be the factor that multiplyin the square of velocity  of joint 1
% we do the same also for joint 2
% the factor 1/2 disappear automatically in this differentation


M(1,1) = simplify(diff(T,dq1,2));
M(2,2)=simplify(diff(T,dq2,2));
M(3,3)=simplify(diff(T,dq3,2));

TempB1=simplify(diff(T,dq1));
M(1,2)= simplify(diff(TempB1,dq2));
M(2,1) = M(1,2);
M(1,3) = simplify(diff(TempB1,dq3));
M(3,1)=M(1,3);

TempB2=simplify(diff(T,dq2));
M(2,3) = simplify(diff(TempB2,dq3));
M(3,2)=M(2,3);
M=simplify(M)
pause


%% parametrization
syms a1 a2 a3 a4 a5 real

Mp=[a1 -a4*cos(q3)*cos(q2) a4*sin(q2)*sin(q3);
     -a4*cos(q3)*cos(q2) a2-a5*sin(q3)^2 0;
     a4*sin(q2)*sin(q3) 0 a3;];
     
disp('*Christoffel matrices*')

%% compute the Christoffel Matrix
q=[q1;q2;q3];

% consider the first column of matrix M, so you have M1
% do the same for M2
M1=Mp(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(Mp,q1))
M2=Mp(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(Mp,q2))
M3=Mp(:,3);
C3=(1/2)*(jacobian(M3,q)+jacobian(M3,q)'-diff(Mp,q3))

pause

disp('***robot centrifugal terms (no Coriolis here!)***')

dq=[dq1;dq2;dq3;];
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c3=dq'*C3*dq;
c=[c1;c2;c3]

pause 

disp('*potential energy of link 1*')

%% vector gravity acceleration
g=[0;0;-g0;];

%% compute the potential energy of link 1
U1=-m1*transpose(g)*pc1

pause

disp('*potential energy of link 2*')

%% compute the potential energy of link 2
U2=-m2*transpose(g)*pc2

pause

%% potential energy of link 3
U3=-m3*transpose(g)*pc3

disp('***robot potential energy (due to gravity)***')

%%  total potential energy
U=simplify(U1+U2+U3)

pause

disp('***robot gravity term***')

%% compute G
G=jacobian(U,q)'

Gp=[a1*g0; a4*g0*cos(q2)*cos(q3);a4*g0*sin(q2)*sin(q3);]
pause

disp('***complete dynamic equations***')

%% show complete dynamic equations
Mp*[ddq1; ddq2; ddq3]+c+Gp==[u1 u2 u3]'

disp('***end***')

% end
