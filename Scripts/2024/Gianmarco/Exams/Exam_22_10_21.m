% @ Exam 22/10/021
% % -------- EXERCISE 1 --------

%% ========= Apply the Moving Frames Algorithm =========
syms Ic1xx Ic1xy Ic1xz Ic1yx Ic1yy Ic1yz Ic1zx Ic1zy Ic1zz Ic2xx Ic2yy Ic2zz real
syms l1 l2 real
syms m1 m2 real
syms q1 q2 q1_dot q2_dot q1_ddot q2_ddot real
syms alpha a d theta
syms g real

n_joints = 2;
q = [q1; q2];
qdots = [q1_dot; q2_dot];
qddots = [q1_ddot; q2_ddot];

% alpha, a_i, d_i, theta_i
DH_table = [ -pi/2,  0,    q1, -pi/2;
            -pi/2,  0,    0,    q2 ;
            ];

% Vector rc containing the values ^ir_{ci} that is the distance
% between the i-th CoM from the i-th reference frame
syms rc1x rc1y rc1z rc2x rc2y rc2z dc2 real
rc1 = [0; 0; 0];
rc2 = [dc2; 0; 0];
rc = [rc1, rc2];

% Vector m containing the masses of the links
m = [m1; m2];

% Baricentric Inertia Matrices
Ic1 = diag([Ic1xx, Ic1yy, Ic1zz]);
Ic2 = diag([Ic2xx, Ic2yy, Ic2zz]);

I = [Ic1, Ic2];

% Moving Frames Algorithm
[w, v, vc, T] = moving_frames_algorithm(n_joints, DH_table, qdots, m, rc, 1, I);

% Extract the kinetic energy T
T_total = simplify(T{1} + T{2});

% Computations for the Inertia Matrix
M = compute_inertia_matrix(T_total, qdots, 3);

% Christoffel Matrix
[C, c_sub] = compute_christoffel_matrix(M, q, qdots, 3);
disp("The Christoffel Matrix is:");
disp(simplify(C));

% Gravity Vector
% In order to retrieve the gravity vector, we need to compute the potential energy
% For the potential energy, we first need to retrieve the heights of the CoMs
h1 = q1;
h2 = q1 + dc2*cos(q2);

% Gravity here is positive, hence g = g
g = g;

U1 = - m1 * g * h1;
U2 = - m2 * g * h2;
U = U1 + U2;

G = transpose(simplify(jacobian(U, [q1, q2])));
G = collect(G, g);

% Complete Dynamic Model | M(q) * ddq + C(q, dq) + G(q) = tau
dynamic_model = simplify(M * transpose([q1_ddot, q2_ddot]) + C + G);

% Linear parametrization of the dynamic model
syms a1 a2 a3 real
av1 = m1 + m2;
av2 = dc2*m2;
av3 = m2*dc2^2 + Ic2yy;
av = [av1, av2, av3];
a = [a1, a2, a3];

%% Inertia Matrix from solutions
% M = [m1 + m2, -m2*dc2*sin(q2);
%     -m2*dc2*sin(q2), m2*dc2^2+Ic2yy];

Ma = subs(M, av, a);
Ca = subs(C, av, a);
Ga = subs(G, av, a);

Y_parametrized = simplify(Ma * transpose([q1_ddot, q2_ddot]) + Ca + Ga);
Y_parametrized = collect(Y_parametrized, sin(q2));
Y = extract_Y_v3(Y_parametrized, a);

disp('Dynamic Model parametrized with a:')
disp(Y)

%% =====================================================
% Now try to solve this exercise as the usual way, meaning computing:
    % - Kinetic Energy
    % - Potential Energy
    % - Compute T
    % - Compute M
    % - Compute C
    % - Compute G
    % - Compute U
    % - Compute the Dynamic Model
