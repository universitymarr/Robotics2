% @ Exam 24/07/08
% -------- EXERCISE 1 --------
syms q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real
syms m1 m2 m3 real
syms dc1 dc2 dc3 real
syms g real
syms Ic2 Ic3 real
syms alpha a d theta
syms l2 real % L2 needed for the DH Table

q = [q1; q2; q3];
qdots = [dq1, dq2, dq3];
qddots = [ddq1, ddq2, ddq3];

% Computations for the Kinetic Energy
% ------------------ T1 ------------------
T1 = 0.5 * m1 * dq1^2;

% ------------------ T2 ------------------
% For computing px and py for q2, we need the DH Table (of a PRR Spatial)
n_joints = 3;

% alpha, a_i, d_i, theta_i
DH_table = [ -pi/2  0    q1    0;
             pi/2  0    l2   q2 ;
              0  dc3    0    q3 ;
            ];

DH = DHmatrix(alpha, a, d, theta);
A = cell(1, n_joints);

% Compute each transformation matrix
for i = 1: n_joints
    A{i} = subs(DH, {alpha, a, d, theta}, DH_table(i, :));
end

T = eye(4);

%disp("-------------------------------------------------------------------")
for i = 1 : n_joints
    % Display print statement
    %disp(['0^A_' num2str(i) '(q_' num2str(i) '):']);
    
    % Perform the simplify operation
    T = simplify(T * A{i});
    
    % Display the simplified result
    %disp(T);
end

% output world-endEff T matrix
%disp("-------------------------------------------------------------------")
%disp('W^T_E:')
%disp(T)

% Remember. T is a matrix 4x4
% Output WE position vector | From 1st to 3rd row of 4th column
p = T(1:3, 4);
%disp('Position vector p:')
%disp(p)

% Rotation Matrix -> From 1st to 3rd row of the first three columns
R = T(1:3, 1:3);

%disp("-------------------------------------------------------------------")

% For computing the linear velocity of the second joint, we need to compute the position vector of the second joint.
% We must, though, compute the position vector of the CoM of the second link
px_joint2 = 0;
py_joint2 = dc2;
pz_joint2 = q1;

% Derive the position vector of the CoM of the second link
dpx_joint2 = diff(px_joint2, q1) * dq1 + diff(px_joint2, q2) * dq2;
dpy_joint2 = diff(py_joint2, q1) * dq1 + diff(py_joint2, q2) * dq2;
dpz_joint2 = diff(pz_joint2, q1) * dq1 + diff(pz_joint2, q2) * dq2;

% ===> Linear velocity of the second joint
velocity_joint2 = [dpx_joint2; dpy_joint2; dpz_joint2];
velocity_joint2 = simplify(velocity_joint2.' * velocity_joint2);

% ===> Angular velocity of the second joint
w2 = dq2;

T2 = simplify(0.5 * m2 * velocity_joint2 + 0.5 * Ic2 * w2^2);

% ------------------ T3 ------------------
% For computing the linear velocity of the third joint, we need to compute the position vector of the third joint.
% We must, though, compute the position vector of the CoM of the third link
px_joint3 = dc3*cos(q2)*cos(q3);
py_joint3 = l2 + dc3*sin(q3);
pz_joint3 = q1 - dc3*cos(q3)*sin(q2);

% Derive the position vector of the CoM of the third link
dpx_joint3 = diff(px_joint3, q1) * dq1 + diff(px_joint3, q2) * dq2 + diff(px_joint3, q3) * dq3;
dpx_joint3 = collect(dpx_joint3, dc3);
dpy_joint3 = diff(py_joint3, q1) * dq1 + diff(py_joint3, q2) * dq2 + diff(py_joint3, q3) * dq3;
dpz_joint3 = diff(pz_joint3, q1) * dq1 + diff(pz_joint3, q2) * dq2 + diff(pz_joint3, q3) * dq3;
dpz_joint3 = collect(dpz_joint3, dc3);

% ===> Linear velocity of the second joint
v_joint3 = simplify(transpose(R) * [dpx_joint3; dpy_joint3; dpz_joint3]);
v_joint3 = simplify(norm(v_joint3.' * v_joint3)^2);

% ===> Angular velocity of the second joint
% To express the angular velocity of the joint 3 we can use the formula: S = R_dot' * R
R_dot = diff(R, q1) * dq1 + diff(R, q2) * dq2 + diff(R, q3) * dq3;
skew_matrix = simplify(R_dot * R');

% Now, for extracting the angular velocity, we need to extract the skew-symmetric matrix terms:
% [-Skew23, skew13, -skew21]
w3_0_3 = [-skew_matrix(2,3); skew_matrix(1,3); -skew_matrix(2,1)];
w3 = simplify(R' * w3_0_3);

T3 = simplify(0.5 * m3 * v_joint3 + 0.5 * transpose(w3) * Ic3 * w3);

% Computation of T
T = simplify(T1 + T2 + T3);

% Computations for the Inertia Matrix
M = compute_inertia_matrix(T, [dq1, dq2, dq3], 3);

% Christoffel Matrix
[C, c_sub] = compute_christoffel_matrix(M, [q1, q2, q3], [dq1, dq2, dq3], 3);

% ------------- MOVING FRAMES ALGORITHM -------------
n_joints = 3;

% alpha, a_i, d_i, theta_i
DH_table = [ -pi/2  0    q1    0;
             pi/2  0    l2   q2 ;
              0  dc3    0    q3 ;
            ];

% Vector m containing the masses of the links
m = [m1; m2; m3];

% Vector rc containing the values ^ir_{ci} that is the distance
% between the i-th CoM from the i-th reference frame
rc = [rc1; rc2; rc3];

% Vector I containing the baricentric inertia matrices of the links
I = [Ic1, Ic2, Ic3];

[w_mf, v_mf, vc_mf, T_mf] = moving_frames_algorithm(n_joints, DH_table, qdots, m, rc, 1, I);
% ---------------------------------------------------

% Gravity Vector
% In order to retrieve the gravity vector, we need to compute the potential energy
% For the potential energy, we first need to retrieve the heights of the CoMs
h1 = q1 - dc1;
h2 = q1;
h3 = q1 - dc3*cos(q3)*sin(q2);

U1 = m1 * g * h1;
U2 = m2 * g * h2;
U3 = m3 * g * h3;
U = U1 + U2 + U3;

G = transpose(simplify(jacobian(U, [q1, q2, q3])));
G = collect(G, g);

% Complete Dynamic Model | M(q) * ddq + C(q, dq) + G(q) = tau
dynamic_model = simplify(M * transpose([ddq1, ddq2, ddq3]) + C + G);

% Linear parametrization of the dynamic model
syms a1 a2 a3 real
av1 = m1 + m2 + m3;
av2 = dc3*m3;
av3 = m3 * dc3^2;
av = [av1, av2, av3];
a = [a1, a2, a3];

Ma = subs(M, av, a);
Ca = subs(C, av, a);
Ga = subs(G, av, a);

Y_parametrized = simplify(Ma * transpose([ddq1, ddq2, ddq3]) + Ca + Ga);
Y = extract_Y_v3(Y_parametrized, a);

disp('Dynamic Model parametrized with a:')
disp(Y)

% -------- EXERCISE 2 --------
% Done by formula 2.4.3 on Formulary

% Define symbolic variables
syms l1 l2 real
syms q1(t) q2(t)
syms q1_dot q2_dot real
syms r_dot_dot real

% Define the joint vector and joint velocity vector as functions of time
q = [q1(t); q2(t)];
q_dot = [diff(q1(t), t); diff(q2(t), t)];

% Task: End-effector position in the x-direction
r = l1*cos(q1(t)) + l2*cos(q1(t) + q2(t));

% Compute the task Jacobian matrix J_r(q) and J_r_dot(q, q_dot)
J_r = jacobian(r, q);
J_r_dot = diff(J_r, t); % time derivative of the Jacobian matrix J_r w.r.t. time

% Compute h(q, q_dot) = J_r_dot * q_dot
h_q_qdot = simplify(J_r_dot * q_dot);

% Define the velocity gain matrix Kd
Kd = diag([2, 2]);

% Compute the pseudoinverse of the task Jacobian matrix
J_r_pseudo = simplify(pinv(J_r));

% Compute the joint accelerations q_dot_dot
q_dot_dot = J_r_pseudo * (r_dot_dot - h_q_qdot) + (eye(2) - J_r_pseudo * J_r) * (-Kd * q_dot);
q_dot_dot = subs(q_dot_dot, [diff(q1(t), t), diff(q2(t), t)], [q1_dot, q2_dot]);

% Substitute the given numerical values
q_dot_dot_subs = subs(q_dot_dot, ...
    [l1, l2, q1(t), q2(t), q1_dot, q2_dot, r_dot_dot], ...
    [1,   1,  pi/4, -pi/2,      1,     -1,         1]);

% Display the result
disp('Joint Accelerations q_dot_dot:')
disp(vpa(q_dot_dot_subs, 4))
