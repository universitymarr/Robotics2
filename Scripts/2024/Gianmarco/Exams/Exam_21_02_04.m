% @ Exam 21/02/04
% -------- EXERCISE 1 --------
syms q1 q2 q3 q4 real
syms dq1 dq2 dq3 dq4 real
syms ddq1 ddq2 ddq3 ddq4 real
syms m1 m2 m3 m4 real
syms I1 I2 I3 I4 real
syms l1 real
syms g
syms dc1 dc2 dc3 dc4 real

q = [q1; q2; q3; q4];
dq = [dq1; dq2; dq3; dq4];
ddq = [ddq1; ddq2; ddq3; ddq4];

num_joints = 4;

% T computations

% -------------------- T1 --------------------
% For computing px and py for q1, we need the DH Table (2R Planar)
% dc_i = Distanza dall'inizio del joint i-1 e il centro di massa del link i
% We use dc_i instead of L, because we are considering the distance from the
% joint i-1 to the center of mass of the link i
px1 = simplify(dc1*cos(q1));
py1 = simplify(dc1*sin(q1));

dpx1 = simplify(time_derivate(px1, q, dq));
dpy1 = simplify(time_derivate(py1, q, dq));

w1 = dq1;                       % Angular velocity of link 1
v1 = simplify(dpx1^2 + dpy1^2); % Linear velocity of link 1 (already squared)

T1 = simplify(0.5 * m1 * v1 + 0.5 * I1 * w1^2);

% -------------------- T2 --------------------
px2 = simplify(l1 * cos(q1));
py2 = simplify(l1 * sin(q1));

dpx2 = simplify(time_derivate(px2, q, dq));
dpy2 = simplify(time_derivate(py2, q, dq));

w2 = dq1 + dq2;                 % Angular velocity of link 2
v2 = simplify(dpx2^2 + dpy2^2); % Linear velocity of link 2 (already squared)

T2 = simplify(0.5 * m2 * v2 + 0.5 * I2 * w2^2);

% -------------------- T3 --------------------
% This is a prismatic joint, so no angular velocity is present
px3 = simplify(l1 * cos(q1) + (q3 - dc3) * cos(q1 + q2));
py3 = simplify(l1 * sin(q1) + (q3 - dc3) * sin(q1 + q2));

dpx3 = simplify(time_derivate(px3, q, dq));
dpy3 = simplify(time_derivate(py3, q, dq));

v3 = simplify(dpx3^2 + dpy3^2); % Linear velocity of link 3 (already squared)

T3 = simplify(0.5 * m3 * v3);

% -------------------- T4 --------------------
px4 = simplify(l1 * cos(q1) + q3 * cos(q1 + q2) + dc4 * cos(q1 + q2 + q4));
py4 = simplify(l1 * sin(q1) + q3 * sin(q1 + q2) + dc4 * sin(q1 + q2 + q4));

dpx4 = simplify(time_derivate(px4, q, dq));
dpy4 = simplify(time_derivate(py4, q, dq));

w4 = simplify(dq1 + dq2 + dq4);           % Angular velocity of link 4
v4 = simplify(dpx4^2 + dpy4^2); % Linear velocity of link 4 (already squared)

T4 = simplify(0.5 * m4 * v4 + 0.5 * I4 * w4^2);

% -------------------------------------------

% Computation of T
T = simplify(T1 + T2 + T3 + T4);

% Inertia Matrix and Christoffel Symbols
M = compute_inertia_matrix(T, dq, 3);
[C, c_sub] = compute_christoffel_matrix(M, q, dq, 3);

% Gravity Vector
% In order to retrieve the gravity vector, we need to compute the potential energy
% For the potential energy, we first need to retrieve the heights of the CoMs
h1 = dc1 * sin(q1);
h2 = l1 * sin(q1);
h3 = l1 * sin(q1) + (q3 - dc3) * sin(q1 + q2);
h4 = l1 * sin(q1) + q3 * sin(q1 + q2) + dc4 * sin(q1 + q2 + q4);

% Potential Energy
U1 = m1 * g * h1;
U2 = m2 * g * h2;
U3 = m3 * g * h3;
U4 = m4 * g * h4;
U = simplify(U1 + U2 + U3 + U4);

% Final Gravity Vector
G = transpose(simplify(jacobian(U, q)));
G = collect(G, g);

% Complete Dynamic Model
dynamic_model = simplify(M * ddq + C + G);

% Linear parametrization of the dynamic model -> TBD

% Equilibrium G = 0
eqn = G == 0;
solutions = solve(eqn, q);
