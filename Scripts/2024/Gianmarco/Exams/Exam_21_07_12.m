% % @ Exam 21/07/12
% % -------- EXERCISE 1 --------
syms q1(t) q2(t) q3(t)
syms L real

% Define the joint vector and joint velocity vector as functions of time
q = [q1(t); q2(t); q3(t)];

% Direct kinematics of 3R robot
p = [L * cos(q1) + L * cos(q1 + q2) + L * cos(q1 + q2 + q3);
     L * sin(q1) + L * sin(q1 + q2) + L * sin(q1 + q2 + q3)];

% The general kinematic control law q_dot is:
J = simplify(jacobian(p, q));
J_dot = diff(J, t);

p_dot = diff(p, t);
p_dot_dot = diff(p_dot, t);

disp('p_dot = ')
disp(p_dot)

disp('p_dot_dot = ')
disp(p_dot_dot)

%% -------- EXERCISE 3 --------
syms q1 q2 q3 real
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms m1 m2 m3 real
syms I1 I2 I3 real
syms l1 real
syms g
syms dc1 dc2 dc3 real

q = [q1; q2; q3];
dq = [dq1; dq2; dq3];
ddq = [ddq1; ddq2; ddq3];

num_joints = 3;

% T computations

% -------------------- T1 --------------------
% For computing px and py for q1, we need the DH Table (2R Planar)
% dc_i = Distanza dall'inizio del joint i-1 e il centro di massa del link i
% We use dc_i instead of L, because we are considering the distance from the
% joint i-1 to the center of mass of the link i
px1 = 0;
py1 = 0;

dpx1 = simplify(time_derivate(px1, q, dq));
dpy1 = simplify(time_derivate(py1, q, dq));

w1 = dq1;                       % Angular velocity of link 1
v1 = simplify(dpx1^2 + dpy1^2); % Linear velocity of link 1 (already squared)

T1 = simplify(0.5 * m1 * v1 + 0.5 * I1 * w1^2);

% -------------------- T2 --------------------
px2 = simplify(q2 - dc2) * cos(q1);
py2 = simplify(q2 - dc2) * sin(q1);

dpx2 = simplify(time_derivate(px2, q, dq));
dpy2 = simplify(time_derivate(py2, q, dq));

w2 = dq1;                 % Angular velocity of link 2
v2 = simplify(dpx2^2 + dpy2^2); % Linear velocity of link 2 (already squared)

T2 = simplify(0.5 * m2 * v2);

% -------------------- T3 --------------------
% This is a prismatic joint, so no angular velocity is present
px3 = simplify(q2 * cos(q1) + dc3*cos(q1 + q3));
py3 = simplify(q2 * sin(q1) + dc3*sin(q1 + q3));

dpx3 = simplify(time_derivate(px3, q, dq));
dpy3 = simplify(time_derivate(py3, q, dq));

w3 = dq1 + dq3;           % Angular velocity of link 3
v3 = simplify(dpx3^2 + dpy3^2); % Linear velocity of link 3 (already squared)

T3 = simplify(0.5 * m3 * v3 + 0.5 * I3 * w3^2);

% -------------------------------------------

% Computation of T
T = simplify(T1 + T2 + T3);

% Inertia Matrix and Christoffel Symbols
M = compute_inertia_matrix(T, dq, 3);
[C, c_sub] = compute_christoffel_matrix(M, q, dq, 3);

% --------------- MOVING FRAMES ALGORITHM ---------------
% alpha, a_i, d_i, theta_i
DH_table = [ -pi/2  0    l1    q1;
              0  l2    0    q2 ;
            ];

% Vector rc containing the values ^ir_{ci} that is the distance
% between the i-th CoM from the i-th reference frame
syms rc1x rc1y rc1z rc2x real
rc1=[rc1x; rc1y; rc1z];
rc2=[rc2x; 0; 0];
rc = [rc1, rc2];

% Vector m containing the masses of the links
m = [m1; m2];

% Baricentric Inertia Matrices
Ic1 = [Ic1xx Ic1xy Ic1xz;
    Ic1xy Ic1yy Ic1yz;
    Ic1xz Ic1yz Ic1zz];
Ic2 = diag([Ic2xx Ic2yy Ic2zz]);

I = [Ic1,Ic2];

% Moving Frames Algorithm
[w_mf, v_mf, vc_mf, T_mf] = moving_frames_algorithm(n_joints, DH_table, qdots, m, rc, [], I);
% -------------------------------------------------------

% Gravity Vector
G = 0;

% Complete Dynamic Model
dynamic_model = simplify(M * ddq + C + G);

% Linear parametrization of the dynamic model
syms a1 a2 a3 a4 a5 a6 a7 a8 a9 real
av1 = dc3 * m3;
av2 = m3 + m2;
av3 = I3 + dc3^2 * m3;
av4 = I1 + I2 + I3 + m3 * dc3;
av5 = dc3;
av6 = dc2;
av7 = I3;
av8 = m3;
av9 = m2;

av = [av1, av2, av3, av4, av5, av6, av7, av8, av9];
a = [a1, a2, a3, a4, a5, a6, a7, a8, a9];

Ma = subs(M, av, a);
Ca = subs(C, av, a);

Ya = simplify(Ma * transpose([ddq1, ddq2, ddq3]) + Ca);
Ysol = extract_Y_v3(Ya, a);

% To verify the correctness of the computation,
% We need to verify that: Ysol * transpose(a) == Ya;