% % @ Exam 22/09/09
% % -------- EXERCISE 1 --------
% syms l1 l2 l3 real
% syms q1(t) q2(t) q3(t)
% syms q1_dot q2_dot q3_dot real
% num_joints = 3;

% % Define the joint vector and joint velocity vector as functions of time
% q = [q1(t); q2(t); q3(t)];

% % Task: End-effector position
% r = [l1 * cos(q1(t)) + l2 * cos(q1(t) + q2(t)) + l3 * cos(q1(t) + q2(t) + q3(t));
%      l1 * sin(q1(t)) + l2 * sin(q1(t) + q2(t)) + l3 * sin(q1(t) + q2(t) + q3(t));];

% % The general kinematic control law q_dot is:
% % q_dot = J# * r_dot + (I - J# * J) * q_dot_0
% % where q_dot_0 = - grad(H(q))
% J = simplify(jacobian(r, q));
% J_pinv = simplify(pinv(J));

% % r_dot
% r_dot = diff(r, t);

% % I
% I = eye(3);

% % H(q)
% q_min = [-(3*pi)/4; -(3*pi)/4; -(3*pi)/4];
% q_max = [(3*pi)/4; (3*pi)/4; (3*pi)/4];
% q_bar = (q_min + q_max) / 2;
% H = sum((q - q_bar).^2 / (q_max - q_min).^2, "all");
% H = 1/(2*num_joints) * H;

% % grad(H)
% gradH = gradient(H, q);

% % q_dot: J# * r_dot + (I - J# * J) * q_dot_0
% % Since we do not move the end-effector, r_dot = 0, hence:
% % q_dot = (I - J# * J) * q_dot_0
% q_dot = simplify((I - J_pinv * J) * -gradH);

% disp('q_dot = (I - J^# * J) * q_dot_0')
% disp(q_dot)

% H_sol = double(subs(H, q, [0, 2*pi/3, 2*pi/3].'));
% grad_H_sol = double(subs(gradH, q, [0, 2*pi/3, 2*pi/3].'));

% disp('H(q) = ')
% disp(H_sol)
% disp('grad(H) = ')
% disp(grad_H_sol)

% -------- EXERCISE 2 --------
% :(

% -------- EXERCISE 3 --------
syms B M g f real
syms theta dtheta ddtheta theta_desired
syms Kp
% Derive the Dynamic model of the system of masses
dynamic_fig_a = (B + M) * ddtheta + (M * g) == f;

lyap_candidate = 1/2 * (B + M) * dtheta^2 + 1/2 * Kp * (theta - theta_desired)^2;
lyap_candidate_dot = diff(lyap_candidate, theta) * dtheta + diff(lyap_candidate, dtheta) * ddtheta;