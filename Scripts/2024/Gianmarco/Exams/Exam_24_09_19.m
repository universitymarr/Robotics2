% EXAM: 19/09/2024
% Student: Scarano Gianmarco - 2047315

% --------- EXERCISE 1 ----------
syms q1 q2 q1_dot q2_dot q1_dot_dot q2_dot_dot real
syms m1 real
syms m2 real
syms dc1 real
syms Ic1 real
syms l1 l2 real
syms g real

n_joints = 2;
q = [q1; q2];
qdots = [q1_dot; q2_dot];
qddots = [q1_dot_dot; q2_dot_dot];

% Computations for the Kinetic Energy
% ------------------ T1 ------------------
% We must compute the position vector of the CoM of the first link
px_joint1 = dc1 * cos(q1);
py_joint1 = dc1 * sin(q1);

% Derive the position vector of the CoM of the first link
dpx_joint1 = diff(px_joint1, q1) * q1_dot;
dpy_joint1 = diff(py_joint1, q1) * q1_dot;

% ===> Linear velocity of the second joint
velocity_joint1 = [dpx_joint1; dpy_joint1];
velocity_joint1 = simplify(velocity_joint1.' * velocity_joint1);

% ===> Angular velocity of the second joint
w1 = q1_dot;

T1 = 0.5 * m1 * velocity_joint1 + 0.5 * w1.' * Ic1 * w1;

% ------------------ T2 ------------------
% Since it is prismatic, we only need to compute the linear velocity of the second joint
% We must compute the position vector of the CoM of the second link
px_joint2 = l1 * cos(q1) + q2 * cos(q1);
py_joint2 = l1 * sin(q1) + q2 * sin(q1);

% Derive the position vector of the CoM of the second link
dpx_joint2 = diff(px_joint2, q1) * q1_dot + diff(px_joint2, q2) * q2_dot;
dpy_joint2 = diff(py_joint2, q1) * q1_dot + diff(py_joint2, q2) * q2_dot;

% ===> Linear velocity of the second joint
velocity_joint2 = [dpx_joint2; dpy_joint2];
velocity_joint2 = simplify(velocity_joint2.' * velocity_joint2);

T2 = 0.5 * m2 * velocity_joint2;

% Computation of T
T = simplify(T1 + T2);

% Computations for the Inertia Matrix
M = compute_inertia_matrix(T, [q1_dot, q2_dot], 3);

% Christoffel Matrix
[C, c_sub] = compute_christoffel_matrix(M, [q1, q2], [q1_dot, q2_dot], 3);

% Gravity Vector
% In order to retrieve the gravity vector, we need to compute the potential energy
% For the potential energy, we first need to retrieve the heights of the CoMs
h1 = dc1 * sin(q1);
h2 = l1 * sin(q1) + q2 * sin(q1);

U1 = m1 * g * h1;
U2 = m2 * g * h2;
U = U1 + U2;

G = transpose(simplify(jacobian(U, [q1, q2])));
G = collect(G, g);

% Complete Dynamic Model | M(q) * ddq + C(q, dq) + G(q) = tau
dynamic_model = simplify(M * transpose([q1_dot_dot, q2_dot_dot]) + C + G);

% Try to expand G and C
G = expand(G);
C = expand(C);

% Linear parametrization of the dynamic model
syms a1 a2 a3 a4 a5 real
av1 = m2;
av2 = m1*dc1;
av3 = l1*m2;
av4 = 2 * l1^2;
av5 = 4 * l1;
av = [av1, av2, av3, av4, av5];
a = [a1, a2, a3, a4, a5];

Ma = subs(M, av, a);
Ca = subs(C, av, a);
Ga = subs(G, av, a);

Y_parametrized = simplify(Ma * transpose([q1_dot_dot, q2_dot_dot]) + Ca + Ga);
Y = extract_Y_v3(Y_parametrized, a);

disp("Inertia Matrix:")
disp(M)

disp("Christoffel Matrix:")
disp(C)

disp("Gravity Vector:")
disp(G)

disp('Dynamic Model:')
disp(dynamic_model)

disp('Dynamic Model parametrized with a:')
disp(Y)

% --------- EXERCISE 2 ----------
syms a b k T real positive
syms t real
syms t_dot real
syms m11 m12 m21 m22 real

M_sym = [m11, m12; m21, m22];
q1_t = a + b * (1 - cos((pi*t)/T));
q2_t = k;

q1_t_dot = time_derivate(q1_t, t, t_dot);
q2_t_dot = time_derivate(q2_t, t, t_dot);

m_multiplied = M_sym * [q1_t_dot; q2_t_dot];
m_multiplied = simplify(m_multiplied);
m_multiplied = subs(m_multiplied, [m11, m21], [Ma(1, 1), Ma(2, 1)]);

% --------- EXERCISE 4 ----------
syms Kd Kp X_eq X_e Y_eq Y_e Fd real positive

eq = Kd*(X_eq - X_e) == Fd;
eq1 = Kp*(Y_eq - Y_e) == 0;

% Isolate now X_eq and Y_eq
X_eq = solve(eq, X_eq);
Y_eq = solve(eq1, Y_eq);
disp('Equilibrium points:')
disp(X_eq)
disp(Y_eq)

% --------- EXERCISE 5 ----------
% Define symbolic variables
syms t a b k pi T real
syms tau1max real

% Define q1(t)
q1_t = a + b * (1 - cos(pi * t / T));
q2_t = k;

q1_dot_t = diff(q1_t, t); % First derivatives (qd_dot)
q1_dot_dot_t = diff(q1_dot_t, t); % Second derivatives (qd_dot_dot)

q2_dot_t = 0; % First derivatives (qd_dot)
q2_dot_dot_t = 0; % Second derivatives (qd_dot_dot)

% Substituting the trajectory into M, C, and G
M_qd = subs(M, [q1, q2], [q1_t, q2_t]);
C_qd = subs(C, [q1, q2, q1_dot, q2_dot], [q1_t, q2_t, q1_dot_t, q2_dot_t]);
G_qd = subs(G, [q1, q2], [q1_t, q2_t]);

% Substitute the second derivatives into the dynamic model
tau_d = M_qd * [q1_dot_dot_t; q2_dot_dot_t] + C_qd * [q1_dot_t; q2_dot_t].' + G_qd;

% Simplify the result
tau_d = simplify(tau_d);

% Extract tau_1 and tau_2
tau_1_t = tau_d(1);
tau_2_t = tau_d(2);
disp('tau_1(t):')
disp(tau_1_t)
disp('tau_2(t):')
disp(tau_2_t)