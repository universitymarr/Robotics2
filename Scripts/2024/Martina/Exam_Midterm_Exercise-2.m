clear all
clc
syms Ixx Iyy Izz q1 q2 q3 q4 g1 g2 qd1 qd2 qd3 qd4 real
syms l1 l2 l3 l4 dc1 dc2 dc3 dc4 m1 m2 m3 m4 vx vy real
syms qdd1 qdd2 qdd3 qdd4 deltaq T l pdx pdy I1 I2 I3 d k tau1 tau2 real
syms rc1x rc1y rc1z rc2x real
syms a1 a2 a3 a4 a5 real;
syms I1xx I1xy I1xz I1yx I1yy I1yz I1zx I1zy I1zz I2xx I2yy I2zz real
digits(4)
%%
DH=[-pi/2,0,l1,q1;
    0,l2,0,q2]
r=[[0;-l1;0],[l2;0;0]]
rc=[[rc1x;rc1y;rc1z],[rc2x;0;0]]
m=[m1;m2]
qd=[qd1;qd2]
I1=[I1xx,I1xy,I1xz;
    I1yx,I1yy,I1yz;
    I1zx,I1zy,I1zz]
I2=[I2xx,0,0;
    0,I2yy,0;
    0,0,I2zz]
I=[I1,I2]
[w,v,vc,T]=moving_frames_algorithm(2,DH,qd,m,r,rc,[],I)
%%
T=T{1}+T{2}
M=simplify(inertia_matrix_from_kinetic_energy(T,qd))
%%

M=[a1+a2*cos(q2)^2+a3*sin(q2)^2,0;
0,a4]
[c,C]=inertia_matrix_to_coriolis(M,[q1,q2],[qd1;qd2])
U1=m1*[-9.81;0;0]'*[rc1x;rc1y;rc1z]
U2=m2*([cos(q2),-sin(q2),0; sin(q2),cos(q2),0;0,0,1]'*[-9.81;0;0])'*[rc2x;0;0]
U=U1+U2
g=[diff(U,q1);diff(U,q2)]
%%
S1=qd'*C{1}
S2=qd'*C{2}
S=[S1;S2]
Mdot=diff(M,q1)*qd1+diff(M,q2)*qd2
skewsymm=Mdot-2*S
g=[0;a5*sin(q2)]
tau=M*[qdd1;qdd2]+c+g
Y=regressor_matrix(tau,[a1;a2;a3;a4;a5])
syms qdd1 qdd2 t real
tau=subs(tau,[q1,q2,qd1,qd2,qdd1,qdd2],[2*t,pi/4,2,0,0,0])