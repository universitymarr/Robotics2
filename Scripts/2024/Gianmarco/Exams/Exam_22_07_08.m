% @ Exam 22/07/08
% -------- EXERCISE 1 --------
% :(

% -------- EXERCISE 2 --------
% :(

% -------- EXERCISE 3 --------
syms q1 q2
syms dq1 dq2
syms qd t T pi
syms a b c d
syms m1 m2 g0
syms ddx ddy
syms ux uy

% General parameters
q = [q1; q2];
dq = [dq1; dq2];
tau = t/T;

% Initial and final conditions
q_in = [1; 0.3];
q_d = [0.6; 0.7];
q_delta = q_d - q_in;
q_delta_x = q_delta(1);
q_delta_y = q_delta(2);

dq_in = [0; 0];
dq_d = [0; 0];

% Dynamic model of the robot
% This must be divided by two equations since we want to retrieve a Dynamic Model
% for each joint (2 prismatic joints)
dynamic_joint_1 = m1 + m2 * ddx == ux;
dynamic_joint_2 = m2 * ddy + m2 * g0 == uy;

% Solve the system by isolating ddx and ddy (we do this manually, hence)
ddx = ux/(m1 + m2);
ddy = (uy/m2) - g0;
ddy_pos = (uy/m2) + g0;

syms ty_sw ty real
q_delta_y = (1/2 * ddy * ty_sw^2) + (1/2 * ddy * ty_sw * (ty - ty_sw));

% We want to isolate ty_sw and ty from the previous equation
ty = sqrt(q_delta_y/ddy) + sqrt(q_delta_y/ddy_pos);
% ty_sw_p = isolate(q_delta_y == 0, ty_sw);
ty_sw_p = ty/2 * ((m2 * g0)/uy);

% Calculations for Tx
tx = 2 * sqrt(abs(q_delta_x/ddx));
tx_sw = tx/2;

% % Polynomial trajectory
% q_pos = q_in + q_delta * (a*tau^3 + b*tau^2 + c*tau + d);
% q_vel= diff(q_pos, t)/T;
% q_acc = diff(q_vel, t)/T;

% cond1 = subs(q_pos, t, 0) == q_in;
% cond2 = subs(q_vel, t, 0) == dq_in;
% cond3 = subs(q_pos, t, T) == q_d;
% cond4 = subs(q_vel, t, T) == dq_d;

% sol = solve([cond1; cond2; cond3; cond4], [a, b, c, d]);
% coef_symbols =  [ a b c d ];
% coef_values = [sol.a sol.b sol.c sol.d];

% q_pos = subs(q_pos, coef_symbols, coef_values );
% q_vel = subs(q_vel, coef_symbols, coef_values);
% q_acc = subs(q_acc, coef_symbols, coef_values);
% q_acc = collect( expand(q_acc), pi);

% syms tau1 tau2 d2 tau_1m tau_2m

% tau = [tau1 tau2];
% tau_m = [tau_1m tau_2m];
% M = [m1 + m2 0 ; 0 m2];

% U1 = 0;
% U2 = m2g0(q2-d2);
% U = U1+U2;
% U = transpose( jacobian(U, q) );

% L = norm(q_delta);
% ddq = inv(M)*( transpose(tau) - U  );
% ddq = subs(ddq, [tau1 tau2], [tau_1m tau_2m])
% ddq_m = subs(ddq, [m1 m2 g0 tau_1m tau_2m], [5 3 9.8 40 40] );
% ddq_m = double(ddq_m)
% T1 = sqrt(abs(4*(q_delta(1)/ddq_m(1))))
% T2 = sqrt(abs(4*q_delta(2)/ddq_m(2)))