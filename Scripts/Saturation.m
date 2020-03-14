%% Saturation for 3R

%% Value of qd
qd=[0; 0; 0;];

%% Bounds
b1 = 2.8;
b2 = 3.6;
b3 = 4;

%% Subs value 
J = simplify(subs(J,{q1,q2,q3},{pi/6,pi/6,pi/6}))
pause
Jd_3R = simplify(subs(Jd_3R,{q1,q2,q3,qd1,qd2,qd3},{pi/6,pi/6,pi/6,qd(1),qd(2),qd(3)}))
pause
rho = rank(J);

%% Compute qdd_ps
qdd_ps= simplify(pinv(J))* (pdd - Jd_3R*qd);

%% Adjust bound limit 
pdd_new = pdd - J(1:2,3)*(-b3);

%% Compute q_new,where it contains it contains the variables
%% that we have not adjusted, but that we need to recalculate
% use pinv because in this way if J is not square but it has full
% row rank, this function compute pseudoinverse
q_new = simplify(pinv(J(1:2,1:2))*pdd_new);

%% Final value qdd and norm
qdd = [-b3; q_new]
norm_qdd = simplify(norm(qdd))
