%% dynamic model of 2R planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 real
syms I1xx I1yy I1zz real
syms I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms d1 d2 l1 l2
syms q1 q2  real
syms dq1 dq2  real 
syms ddq1 ddq2 real 
syms u1 u2 g0 real

disp('**** dynamic model of 2R planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

%% initialization
v0=0;
w0=0;
%% Rotation Matrix R
R1 = [cos(q1) -sin(q1) 0;
      sin(q1) cos(q1) 0;
      0 0 1]

R2 = [cos(q2) -sin(q2) 0;
      sin(q2) cos(q2) 0;
      0 0 1]

%% position of center of mass of each link  
rc1 = [-l1+d1; 0; 0;];
rc2 = [-l2+d2; 0; 0;];

%% position
pc1=[d1*cos(q1);d1*sin(q1); 0];
pc2=[l1*cos(q1) + d2*cos(q1+q2); l1*sin(q1)+d2*sin(q1+q2);0];


disp(' ')
disp('*START RECURSION OF ALGORITHM*')
disp(' ')
disp('*angular velocity of link 1 with respect to frame 1*')

%% compute angular velocity w of link 1
w1 = transpose(R1)*(w0+dq1*[0;0;1]) 

disp(' ')
disp('*linear velocity of link 1 with respect to frame 1*')

%% compute linear velocity v of link 1
v1 = transpose(R1)*(v0+cross(w1,[l1*cos(q1);l1*sin(q1);0]))

disp(' ')
disp('*linear velocity of center od mass of link 1 with respect to frame 1*')

%% compute velocity of center of mass of link 1
vc1 = v1 + cross(w1,rc1)

disp(' ')
disp('*kinetic energy of link 1*')

%% compute kinetic energy of joint 1
T1=(1/2)*(I1zz + m1*d1^2)*dq1^2

pause

disp(' ')
disp('*angular velocity of link 2 with respect to frame 2*')

%% compute angular velocity w of link 1
w2 = transpose(R2)*(w1+dq2*[0;0;1]) 

disp(' ')
disp('*linear velocity of link 2 with respect to frame 2*')

%% compute linear velocity v of link 1
v2 = transpose(R2)*(v1+cross(w2,[l2*cos(q2);l2*sin(q2);0]))

disp(' ')
disp('*linear velocity of center od mass of link 1 with respect to frame 1*')

%% compute velocity of center of mass of link 1
vc2 = v2 + cross(w2,rc2)

disp('*kinetic energy of link 2 *')

%% kinetic energy of joint 2
T2=1/2*m2*(l1^2*dq1^2+d1^2*(dq1+dq2)^2+2*l1*d2*cos(q2)*dq1*(dq1+dq2))+ 1/2*I2zz*(dq1+dq2)^2;
pause

disp('***robot kinetic energy***')

%% total kinetic energy
T=T1+T2;

disp('*simplifying*')

T=simplify(T1+T2)
pause

% collect in base of term that you pass it in this case, collect terms that has dq1^2 and 
% do the same for dq2^2 because you need them to compute M
T=collect(T,dq1^2)
T=collect(T,dq2^2)

pause

disp('***robot inertia matrix***')

%% compute robot matrix M, inertia matrix

% we extract the element M(1,1) from the total kinetic energy by doing twice the derivative
% respect qd1 because this will be the factor that multiplyin the square of velocity  of joint 1
% we do the same also for joint 2
% the factor 1/2 disappear automatically in this differentation

a1 = I1zz + m1*d1^2 +I2zz+m2*d2^2+m2*l1^2;
a2 = m2*l1*d2;
a3 = I2zz+m2*d2^2;

M(1,1) = a1+2*a2*cos(q2)
M(1,2)= a3+a2*cos(q2);
M(2,1) = a3+a2*cos(q2);
M(2,2)=a3;
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

%% vector gravity acceleration
g=[0;-g0;0];

%% compute the potential energy of link 1
U1=-m1*transpose(g)*pc1

pause

disp('*potential energy of link 2*')

%% compute the potential energy of link 2
U2=-m2*transpose(g)*pc2

pause

disp('***robot potential energy (due to gravity)***')

%%  total potential energy
U=simplify(U1+U2)

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
