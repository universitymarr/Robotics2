%% dynamic model of 2R motor planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 real
syms mm1 mm2 real
syms L1 L2 real
syms d1 d2 real
syms I1xx I1yy I1zz real  %symbolic variables explicitly defined as real
syms Im1xx Im1yy Im1zz real  %symbolic variables explicitly defined as real
syms I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms Im2xx Im2yy Im2zz real  %symbolic variables explicitly defined as real
syms t1 t2 real
syms tm1 tm2 real
syms dt1 dt2 real
syms dtm1 dtm2 real
syms ddt1 ddt2 real
syms ddtm1 ddtm2 real
syms u1 u2  u3 u4 g0 real

disp('**** dynamic model of 2R planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

disp(' ')
disp('*kinetic energy of link 1 - linear part*')

%% compute linear part of kinetic energy of joint 1
pc1 =[d1*cos(t1) d1*sin(t1) 0]';
vc1 = diff(pc1,t1)*dt1;
Tl1 = (1/2)*m1*vc1'*vc1;

disp('*kinetic energy of link 1 - angular part*')

%% compute the angular part of kinetic energy of joint 2
%om1 = [0 dq1 0]';
om1 = [0 0 dt1]';
Ta1 = (1/2)*om1'*diag([I1xx I1yy I1zz])*om1;

T1= simplify(Tl1+Ta1)
pause

disp('*kinetic energy of link 2 - linear part*')

%% compute the linear part of kinetic energy of joint 2
%pc2 = [d2*cos(q1)*cos(q2) d2*cos(q2)*sin(q1) L1+d2*sin(q2)]';
pc2 = [L1*cos(t1)+d2*cos(t2) L1*sin(t1)+d2*sin(t2) 0]';
vc2 = simplify(diff(pc2,t1)*dt1+diff(pc2,t2)*dt2);
T2l = simplify((1/2)*m2*vc2'*vc2);

disp('*kinetic energy of link 2 - angular part*')

%% compute the angular part of kinetic energy of joint 2
%om2=[dq1*sin(q2) dq1*cos(q2) dq2]';
om2=[0 0 dt2]';
T2a=simplify((1/2)*om2'*diag([I2xx I2yy I2zz])*om2);

%% kinetic energy of joint 2
T2=simplify(T2l+T2a)
pause

disp('***robot motor energy***')
omm1=[0 0 dtm1]';
Tm1 = simplify((1/2)*omm1'*diag([Im1xx Im1yy Im1zz])*omm1)
pause

omm2=[0 0 dtm2+dt1]';
Tm2 = simplify((1/2)*omm2'*diag([Im2xx Im2yy Im2zz])*omm2) + (1/2)*mm2*L1^2*dt1^2
pause

disp('***robot kinetic energy***')
    
%% 
%% total kinetic energy
T=T1+T2+Tm1+Tm2;

disp('*simplifying*')
    
T=simplify(T);

% collect in base of term that you pass it in this case, collect terms that has dq1^2 and 
% do the same for dq2^2 because you need them to compute M
T=collect(T,dt1^2);
T=collect(T,dt2^2);
T=collect(T,dtm1^2);
T=collect(T,dtm2^2)
pause

disp('***robot inertia matrix***')

%% compute robot matrix M, inertia matrix

% we extract the element M(1,1) from the total kinetic energy by doing twice the derivative
% respect qd1 because this will be the factor that multiplyin the square of velocity  of joint 1
% we do the same also for joint 2
% the factor 1/2 disappear automatically in this differentation
M(1,1)=diff(T,dt1,2);
M(2,2)=diff(T,dt2,2);
M(3,3)= diff(T,dtm1,2);
M(4,4)=diff(T,dtm2,2);

TempB1=diff(T,dt1);
M(1,2)=diff(TempB1,dt2);
M(2,1)=M(1,2);

M(1,3)=diff(TempB1,dtm1);
M(3,1)=M(1,3);

M(1,4)=diff(TempB1,dtm2);
M(4,1)=M(1,4);

TempB2=diff(T,dt2);
M(2,3) = diff(TempB2,dtm1);
M(3,2)=M(2,3);

M(2,4)=diff(TempB2,dtm2);
M(4,2)=M(2,4);

TempB3=diff(T,dtm1);
M(3,4)=diff(TempB3,dtm2);
M(4,3)=M(3,4);

M=simplify(M)

pause

%% parametrization 
syms a1 a2 a3 a4 real
% M=[ a1 a3*cos(q1-q2);
%     a3*cos(q1-q2) a2;]
disp('*Christoffel matrices*')

%% compute the Christoffel Matrix
q=[t1;t2;tm1;tm2];

% consider the first column of matrix M, so you have M1
% do the same for M2
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,t1))
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,t2))
M3=M(:,1);
C3=(1/2)*(jacobian(M3,q)+jacobian(M3,q)'-diff(M,tm1))
M4=M(:,2);
C4=(1/2)*(jacobian(M4,q)+jacobian(M4,q)'-diff(M,tm2))
pause

disp('***robot centrifugal terms (no Coriolis here!)***')

dq=[dt1;dt2;dtm1;dtm2];
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
U1=m1*g'*[0;d1;0]

pause

disp('*potential energy of link 2*')

%% compute the potential energy of link 2
U2=m2*g'*[0;d2*sin(q2);0]
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
