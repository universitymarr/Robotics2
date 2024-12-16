%% Apply the Moving Frames Algorithm
syms Ic1xx Ic1xy Ic1xz Ic1yx Ic1yy Ic1yz Ic1zx Ic1zy Ic1zz Ic2xx Ic2yy Ic2zz real
syms l1 l2 real
syms m1 m2 real
syms q1 q2 q1_dot q2_dot real
syms alpha a d theta

n_joints = 2;
q = [q1; q2];
qdots = [q1_dot; q2_dot];

% alpha, a_i, d_i, theta_i
DH_table = [ -pi/2,  0,    l1,    q1;
              0,  l2,    0,    q2 ;
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
Ic2 = diag([Ic2xx, Ic2yy, Ic2zz]);

I = [Ic1, Ic2];

% Moving Frames Algorithm
[w, v, vc, T] = moving_frames_algorithm(n_joints, DH_table, qdots, m, rc, [], I);

% Extract the kinetic energy T
T = simplify(T{1} + T{2});

% Computations for the Inertia Matrix
M = compute_inertia_matrix(T, qdots, 3);

% EXTRA: Substitutions for simplifications
M = subs(M, sin(q2)^2, (1-cos(q2)^2));
M = simplify(M);
M = collect(M, cos(q2));
M = simplify(M);

disp("The Inertia Matrix is:");
disp(M);

% Christoffel Matrix
[C, c_sub] = compute_christoffel_matrix(M, q, qdots, 3);
disp("The Christoffel Matrix is:");
disp(simplify(C));

% Linear parametrization
% In this case we already know the values from the solutions
% but we should have parametrized the values manually
syms a1 a2 a3 real
M = [a1+a2*cos(q2)^2, 0;
                   0, a3];
[c,C] = compute_christoffel_matrix(M, q, qdots, 3);