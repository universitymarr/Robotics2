% @ Exam 21/07/12
% -------- EXERCISE 1 --------
syms a2 a3 real positive
syms rc1x rc1y rc1z rc2x real
syms m1 m2 m3 real positive
syms Ic1xx Ic1yy Ic1zz Ic2xx Ic2yy Ic2zz Ic3xx Ic3yy Ic3zz real positive
syms dc1 dc2 dc3 real positive
syms theta_1 theta_2 theta_3 real
syms theta_1_dot theta_2_dot theta_3_dot real
syms theta_1_ddot theta_2_ddot theta_3_ddot real

% --------------- MOVING FRAMES ALGORITHM ---------------
n_joints = 3;

theta = [theta_1; theta_2; theta_3];
thetadots = [theta_1_dot; theta_2_dot; theta_3_dot];
thetaddots = [theta_1_ddot; theta_2_ddot; theta_3_ddot];

% alpha, a_i, d_i, theta_i
DH_table = [ pi/2  0   0    theta_1;
              0  a2    0    theta_2;
              0  a3    0    theta_3;
            ];

% Vector m containing the masses of the links
m = [m1; m2; m3];

% Vector rc containing the values ^ir_{ci} that is the distance
% between the i-th CoM from the i-th reference frame
rc1 = [dc1; 0; 0];
rc2 = [dc2 - a2; 0; 0];
rc3 = [dc3 - a3; 0; 0];

rc = [rc1, rc2, rc3];

% Vector I containing the inertia matrices of the links
Ic1 = diag([Ic1xx Ic1yy Ic1zz]);
Ic2 = diag([Ic2xx Ic2yy Ic2zz]);
Ic3 = diag([Ic3xx Ic3yy Ic3zz]);

I = [Ic1, Ic2, Ic3];

[w_mf, v_mf, vc_mf, T_mf] = moving_frames_algorithm(n_joints, DH_table, thetadots, m, rc, [], I);

% Extract the T matrix
T = simplify(T_mf{1} + T_mf{2} + T_mf{3});

% Inertia Matrix and Christoffel Symbols
M = compute_inertia_matrix(T, thetadots, 3);
[C, c_sub] = compute_christoffel_matrix(M, theta, thetadots, 3);

% -------- EXERCISE 2 --------
