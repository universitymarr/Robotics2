clear all;
clc; 

syms q1 q2 q3 l1 l2 l3 m1 m2 m3 dc1 dc2 dc3 real 
digits(4);

g0 = 9.81; 
g = [g0; 0; 0];
rc1 = [-l1+dc1;0;0];
rc2 = [-l2+dc2;0;0];
rc3 = [-l3+dc3;0;0];

DH = [0 l1 0 q1; 0 l2 0 q2; 0 l3 0 q3] 

zero_A1 = DHMatrix(DH(1,:));
zero_R1 = zero_A1(1:3,1:3);
zero_A2 = DHMatrix(DH(1:2,:));
zero_R2 = zero_A2(1:3,1:3);
zero_A3 = DHMatrix(DH);
zero_R3 = zero_A3(1:3,1:3);

one_A2 = DHMatrix(DH(2,:));
one_R2 = one_A2(1:3,1:3);
two_A3 = DHMatrix(DH(3,:));
two_R3 = two_A3(1:3,1:3);

% Forces f3, f2, f1
f3 = vpa(simplify(-m3*zero_R3.'*g)) 
f2 = vpa(simplify(two_R3*f3 -m2*zero_R2.'*g))
f1 = vpa(simplify(one_R2*f2 -m1*zero_R1.'*g))

% Tau3 
three_r23 = [l3;0;0];
tau3 = vpa(simplify(cross(-f3,(three_r23+rc3))))

% Tau2
two_r12 = [l2;0;0];
first_tau2 = vpa(simplify(two_R3*tau3));
second_tau2 = vpa(simplify(cross(two_R3*f3,rc2)));
third_tau2 = vpa(simplify(cross(-f2,(two_r12+rc2))));
tau2 = vpa(simplify(first_tau2+second_tau2+third_tau2))

% Tau1
one_r01 = [l1;0;0];
first_tau1 = vpa(simplify(one_R2*tau2));
second_tau1 = vpa(simplify(cross(one_R2*f2,rc1)));
third_tau1 = vpa(simplify(-cross(f1,(one_r01+rc1))));
tau1 = vpa(simplify(first_tau1 + second_tau1 + third_tau1))

% u3 
three_z2 = vpa(simplify(two_R3.'*[0;0;1]));
u3 = vpa(simplify(tau3.'*three_z2))

% u2 
two_z1 = vpa(simplify(one_R2.'*[0;0;1]));
u2 = vpa(simplify(tau2.'*two_z1))

% u1
one_z0 = vpa(simplify(zero_R1.'*[0;0;1]));
u1 = vpa(simplify(tau1.'*one_z0))