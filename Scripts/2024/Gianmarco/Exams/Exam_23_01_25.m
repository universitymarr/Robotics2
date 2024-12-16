% @ Exam 23/01/25
% -------- EXERCISE 1.a --------
syms q1 q2 real
syms dq1 dq2 real
syms ddq1 ddq2 real
syms m1 m2 I1 I2 real
syms l1 g real
syms dc1 real
syms fv1 fv2 real

q = [q1; q2];
dq = [dq1; dq2];
ddq = [ddq1; ddq2];

num_joints = 2;

% T computations

% -------------------- T1 --------------------
% For computing px and py for q1, we need the DH Table (2R Planar)
% dc_i = Distanza dall'inizio del joint i-1 e il centro di massa del link i
% We use dc_i instead of L, because we are considering the distance from the
% joint i-1 to the center of mass of the link i
px1 = dc1 * cos(q1);
py1 = dc1 * sin(q1);

dpx1 = time_derivate(px1, q, dq);
dpy1 = time_derivate(py1, q, dq);

w1 = dq1; % Angular velocity of link 1
v1 = simplify(dpx1^2 + dpy1^2); % Linear velocity of link 1

T1 = 0.5 * m1 * v1 + 0.5 * I1 * w1^2;

% -------------------- T2 --------------------
% For computing px and py for q2, we need the DH Table (2R Planar)
px2 = l1 * cos(q1);
py2 = l1 * sin(q1);

dpx2 = time_derivate(px2, q, dq);
dpy2 = time_derivate(py2, q, dq);

w2 = dq1 + dq2; % Angular velocity of link 2
v2 = simplify(dpx2^2 + dpy2^2); % Linear velocity of link 2

T2 = 0.5 * m2 * v2 + 0.5 * I2 * w2^2;

% -------------------------------------------

% Computation of T
T = simplify(T1 + T2);

% Inertia Matrix and Christoffel Symbols
M = compute_inertia_matrix(T, dq, 3);
[C, c_sub] = compute_christoffel_matrix(M, q, dq, 3);

% Gravity Vector
% In order to retrieve the gravity vector, we need to compute the potential energy
% For the potential energy, we first need to retrieve the heights of the CoMs
h1 = dc1 * cos(q1);
h2 = l1 * cos(q1);

U1 = -m1 * g * h1;
U2 = -m2 * g * h2;
U = U1 + U2;

G = transpose(simplify(jacobian(U, q)));
G = collect(G, g);

% Friction Vector: -F * q_dot (For a single joint is: -Fv_i * dq_i)
Fv = [fv1; fv2];
F = Fv .* dq;

% Complete Dynamic Model
dynamic_model = simplify(M * ddq + C + G + F);
disp('Dynamic Model:')
disp(dynamic_model)

% Linear parametrization of the dynamic model
syms a1 a2 a3 a4 a5 real
av1 = m1*dc1^2 + m2*l1^2 + I1 + I2;
av2 = dc1*m1 + l1*m2;
av3 = I2;
av4 = fv1;
av5 = fv2;
av = [av1, av2, av3, av4, av5];
a = [a1, a2, a3, a4, a5];

Ma = subs(M, av, a);
Ca = subs(C, av, a);
Ga = subs(G, av, a);

Ya = simplify(Ma * transpose([ddq1, ddq2]) + Ca + Ga);
Ysol = extract_Y_v3(Ya, a);

% -------- EXERCISE 1.b --------
% As for the steps, we need to find G, but we already computed it in the previous step

% We need to find now the Jacobian of G
G_q = simplify(jacobian(G, q));

% Find the maximum eigenvalue of G_q
eig_G_q = simplify(eig(G_q));
max_eig = simplify(max(eig_G_q(2)));
disp('Max eigenvalue of G_q:')
disp(max_eig)

alpha = g * (dc1*m1 + l1*m2);
disp("Lower bound:")
disp(alpha)

% Kp must be >= alpha, but we do not choose it explicitly
% Kd must be 0 instead, since we have viscous friction
% In order to find a proper Kp:
% ||G(q) - G(qd)|| <= Kp * ||q - qd||
% Then solve the inequality for Kp

% -------- EXERCISE 1.c on paper --------