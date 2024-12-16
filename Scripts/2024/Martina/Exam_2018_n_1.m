clear all
clc
syms Ixx Iyy Izz q1 q2 q3 q4 g1 g2 qd1 qd2 qd3 qd4 l1 l2 l3 l4 dc1 dc2 dc3 dc4 m1 m2 m3 m4 real
digits(6)
%%
r=[[0;q1;0],[0;q2;0],[l3;0;0],[l4;0;0]]
rc=[[0;-dc1;0],[0;-dc2;0],[-l3+dc3;0;0],[-l4+dc4;0;0]]
DH=[pi/2,0,q1,pi/2;pi/2,0,q2,pi/2;
0,l3,0,q3;
0,l4,0,q4]
[w,v,vc,T]=moving_frames_algorithm(4,DH,[qd1;qd2;qd3;qd4],[m1;m2;m3;m4],r,rc,[1,2])
T_tot=simplify(T{1}+T{2}+T{3}+T{4})
T=(m4*((dc4*qd3 + dc4*qd4 + qd2*cos(q3 + q4) - qd1*sin(q3 + q4) + l3*qd3*cos(q4))^2 + (sin(q4)*(l3*qd3 + qd2*cos(q3) - qd1*sin(q3)) + cos(q4)*(qd1*cos(q3) + qd2*sin(q3)))^2))/2 + (Izz*qd3^2)/2 + (m3*((qd1*cos(q3) + qd2*sin(q3))^2 + (dc3*qd3 + qd2*cos(q3) - qd1*sin(q3))^2))/2 + (m1*qd1^2)/2 + (m2*(qd1^2 + qd2^2))/2 + Izz*(qd3/2 + qd4/2)*(qd3 + qd4)
M=simplify(jacobian(jacobian(T,[qd1,qd2,qd3,qd4]),[qd1,qd2,qd3,qd4]))