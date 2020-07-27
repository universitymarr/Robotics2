
clear all
close all
clc

%% define symbolic variables
syms m1 m2 m3 m4 real
syms d1 d2 d3 d4 real
syms L1 L2 L3 L4 L h real
syms I1xx I1yy I1zz real  %symbolic variables explicitly defined as real
syms I2xx I2yy I2zz real  %symbolic variables explicitly defined as real
syms I3xx I3yy I3zz real  %symbolic variables explicitly defined as real
syms I4xx I4yy I4zz real  %symbolic variables explicitly defined as real
syms q1 q2 q3 q4 real
syms dq1 dq2 dq3 dq4 real 
syms ddq1 ddq2 ddq3 ddq4 real 
syms u1 u2 u3 u4 g0 real

%% 1
pc1 = [q1; 0;0];
vc1=diff(pc1,q1)*dq1;
T1 = 1/2*m1*vc1'*vc1
pause


%% 2
pc2 = [q1+d2*cos(q2); h+d2*sin(q2);0];
vc2=diff(pc2,q1)*dq1+diff(pc2,q2)*dq2;
Tl2 = 1/2*m2*vc2'*vc2;

om2 = [0 0 dq2]';
Ta2 = 1/2*om2'*diag([I2xx I2yy I2zz])*om2;

T2=Tl2+Ta2
pause


%% 3
pc3 = [q1+L*cos(q2)+d3*cos(q3); h+L*sin(q2)+d3*sin(q3);0];
vc3=diff(pc3,q1)*dq1+diff(pc3,q2)*dq2++diff(pc3,q3)*dq3;
Tl3 = 1/2*m3*vc3'*vc3;

om3 = [0 0 dq3]';
Ta3 = 1/2*om3'*diag([I3xx I3yy I3zz])*om3;

T3=Tl3+Ta3
pause


%% 4

pc4 = [q1+L*cos(q2)+L*cos(q3)+d4*cos(q4); h+L*sin(q2)+L*sin(q3)+d4*sin(q4);0];
vc4=diff(pc4,q1)*dq1+diff(pc4,q2)*dq2+diff(pc4,q3)*dq3+diff(pc4,q4)*dq4;
Tl4 = 1/2*m4*vc4'*vc4

om4 = [0 0 dq4]';
Ta4 = 1/2*om4'*diag([I4xx I4yy I4zz])*om4;

T4= simplify(Tl4+Ta4)
pause

T=simplify(T1+T2+T3+T4)
pause

T=collect(T,dq1^2);
T=collect(T,dq2^2);
T=collect(T,dq3^2);
T=collect(T,dq4^2);


%% inertia

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
