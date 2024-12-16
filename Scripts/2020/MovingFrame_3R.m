%% dynamic model of 3R planar robot under gravity
%% using a Lagrangian formulation in symbolic form

clear all
close all
clc

%% define symbolic variables
syms m1 m2 m3 d 
syms I1xx I1yy I1zz real
syms I2xx I2yy I2zz real
syms I3xx I3yy I3zz real
syms d1 d2 L1 L2 L3 real
syms q1 q2 q3 real
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms u1 u2 u3 g0 real

disp('**** dynamic model of 2R planar robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

%% initialization
w0=[0;0;0];
v0=[0;0;0];
%% R matrix

R1 = [cos(q1) 0  sin(q1);
      sin(q1) 0 -cos(q1);
        0 1 0];

R2 = [cos(q2) -sin(q2) 0;
    
    sin(q2) cos(q2) 0;
    0 0 1;];

R3 = [cos(q3) -sin(q3) 0;
    sin(q3) cos(q3) 0;
    0 0 1;];

%% position of center of mass of each link
syms A F C D E real

rc1 = [A;-F;0;];
rc2 = [-C;0;0];
rc3 = [-D;0;E];

%% positions of link
pc1 = R1'*[0 0 L1]'
pc2 = simplify(R2'*[L2*cos(q2) L2*sin(q2) 0]')
pc3 = simplify(R3'*[L3*cos(q3) L3*sin(q3) 0]')
 
pause
disp('*START RECURSION OF ALGORITHM*')
disp(' ')

%% compute angular velocity w of link 1
w1 = transpose(R1)*(w0+dq1*[0;0;1]);

%% compute linear velocity v of link 1
v1 = transpose(R1)*v0 + cross(w1, pc1);

%% compute velocity of center of mass of link 1
vc1 = v1 + cross(w1,rc1)

disp(' ')
disp('*kinetic energy of link 1*')

%% compute kinetic energy of joint 1
Tl1=(1/2)*m1*vc1'*vc1;
Ta1 = 1/2*w1'*diag([I1xx I1yy I1zz])*w1;

T1 = simplify(Tl1+Ta1)
pause

disp('*angular velocity of link 1 with respect to frame 2*')

%% compute angular velocity w of link 2
w2 = transpose(R2)*(w1+dq2*[0;0;1]); 

%% compute linear velocity v of link 2
v2 = transpose(R2)*v1 +cross(w2, pc2);

%% compute velocity of center of mass of link 2
vc2 = v2 + cross(w2,rc2)
pause
disp('*kinetic energy of link 2*')

%% compute kinetic energy of joint 2
Tl2=(1/2)*m2*vc2'*vc2;
Ta2 = 1/2*w2'*diag([I2xx I2yy I2zz])*w2;

T2=simplify(Tl2+Ta2)
pause

disp('*angular velocity of link 1 with respect to frame 3*')

%% compute angular velocity w of link 3
w3 = transpose(R3)*(w2+dq3*[0;0;1]); 

%% compute linear velocity v of link 2
v3 = transpose(R3)*v2+cross(w3, pc3);

%% compute velocity of center of mass of link 2
vc3 = v3 + cross(w3,rc3)
pause
disp('*kinetic energy of link 3*')

%% compute kinetic energy of joint 2
Tl3=(1/2)*m3*vc3'*vc3;
Ta3 = 1/2*w3'*diag([I3xx I3yy I3zz])*w3;

T3=simplify(Tl3+Ta3)
pause

T=simplify(T1+T2+T3)

T=collect(T,dq1^2);
T=collect(T,dq2^2);
T=collect(T,dq3^2);


%% compute M

M(1,1)=diff(T,dq1,2);
M(2,2)=diff(T,dq2,2);
M(3,3)=diff(T,dq3,2);

TempB12=diff(T,dq1);
M(1,2)=diff(TempB12,dq2);
M(2,1)=M(1,2);

M(1,3)=diff(TempB12,dq3);
M(3,1)=M(1,3);

TempB2=diff(T,dq2);
M(2,3)=diff(TempB2,dq3);
M(3,2)=M(2,3);
M=simplify(M)