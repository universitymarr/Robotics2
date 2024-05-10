% dynamic model of PR planar robot under gravity
% Lagrangian formulation in symbolic form
% A. De Luca
% 23 March 2010

clear all
clc

syms m1 m2 d I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms q1 q2 dq1 dq2 g0 real

disp('**** dynamic model of PR planar robot in the vertical plane ****')
disp(' ')

pause

disp('*kinetic energy of link 1*')

T1=(1/2)*m1*dq1^2

pause

disp('*kinetic energy of link 2 - linear part*')

pc2=[q1+d*cos(q2) d*sin(q2) 0]'
vc2=diff(pc2,q1)*dq1+diff(pc2,q2)*dq2
T2c= (1/2)*m2*vc2'*vc2

pause

disp('*kinetic energy of link 2 - angular part*')

om2=[0 0 dq2]'
T2a=(1/2)*om2'*diag([I2xx I2yy I2zz])*om2

T2=T2c+T2a;
pause

disp('***robot kinetic energy***')

T=T1+T2

pause

disp('*simplifying*')

T=simplify(T1+T2)
T=collect(T,dq1^2)
T=collect(T,dq2^2)

pause

disp('***robot inertia matrix***')

B(1,1)=diff(T,dq1,2);
B(2,2)=diff(T,dq2,2);
TempB12=diff(T,dq1);
B(1,2)=diff(TempB12,dq2);
B(2,1)=B(1,2);
B=simplify(B)

pause

disp('*Christoffel matrices*')

q=[q1;q2];
b1=B(:,1);
C1=(1/2)*(jacobian(b1,q)+jacobian(b1,q)'-diff(B,q1))
b2=B(:,2);
C2=(1/2)*(jacobian(b2,q)+jacobian(b2,q)'-diff(B,q2))

pause

disp('***robot Coriolis and centrifugal term***')

dq=[dq1;dq2];
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c=[c1;c2]

pause 

disp('*potential energy of link 1*')

U1=0

pause

disp('*potential energy of link 2*')

g=[0;-g0;0]
    
U2=-m2*g'*pc2

pause

disp('***robot potential energy (due to gravity)***')

U=simplify(U1+U2)

pause

disp('***robot gravity term***')

G=jacobian(U,q)'

disp('***end***')

% end
