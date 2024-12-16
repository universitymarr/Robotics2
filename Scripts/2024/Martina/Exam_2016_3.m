clear all
clc
syms Ixx Iyy Izz q1 q2 q3 q4 g1 g2 qd1 qd2 qd3 qd4 real
syms l1 l2 l3 l4 dc1 dc2 dc3 dc4 m1 m2 m3 m4 vx vy real
syms a1 a2 a3 a4 a5 a6 th thd m I rc1x rc2x rc1y rc2y I3 I4 real
syms qdd1 qdd2 qdd3 qdd4 deltaq T l pdx pdy I1 I2 I3 d k tau1 tau2 real
syms Iyy1 Ixx2 Iyy2 Izz2 rdot real
digits(6)
%%
p=[cos(q1+q2+q3)+cos(q1+q2)+cos(q1);
    sin(q1+q2+q3)+sin(q1+q2)+sin(q1)]
j=jacobian(p,[q1,q2,q3])
j=subs(j,[q1,q2,q3],[pi/6,pi/6,pi/6])
qd=[0;0;0]
jdot=diff(j,q1)*qd1+diff(j,q2)*qd2+diff(j,q3)*qd3
pddot=[4;2]
qdot=simplify(pinv(j)*(pddot-jdot*qd))
q=vpa(subs(qdot,[q1,q2,q3],[pi/6,pi/6,pi/6]))
pddotnew=pddot-j(:,3)*-4
jnew=j(:,1:2)
qddnew=vpa(subs(pinv(jnew)*(pddotnew-jdot(:,1:2)*qd(1:2)),[q1,q2,q3],[pi/6,pi/6,pi/6]))