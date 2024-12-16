syms q1 q2 q3
syms a2 a3
syms theta1 theta2 theta3
syms dc1

% alpha, a_i, d_i, theta_i
syms alpha a d theta
DH_table = [   pi/2 0 0 theta1;
                0 a2 0 theta2;
                0 a3 0 theta3;
            ];

DH = DHmatrix(alpha, a, d, theta);

A = cell(1,n_joints);

% Compute each transformation matrix
% -- Remember:
%   A{1} = 0^A_1
%   A{2} = 1^A_2
%   A{3} = 2^A_3
%   A{4} = 3^A_4
%   ...
%   A{i} = i-1^A_i
for i = 1: n_joints
    A{i} = subs(DH, {alpha, a, d, theta}, DH_table(i, :));
end

T = eye(4);

disp("-------------------------------------------------------------------")
for i = 1 : n_joints
    % Display print statement
    disp(['0^A_' num2str(i) '(q_' num2str(i) '):']);
    
    % Perform the simplify operation
    T = simplify(T * A{i});
    
    % Display the simplified result
    disp(T);
end

% v = alpha, a d, theta
A01 = A{1};
A12 = A{2};
A23 = A{3};

R01 = A01(1:3, 1:3);
R12 = A12(1:3, 1:3);
R23 = A23(1:3, 1:3);

t1 = A01(1:3, 4);
t2 = A12(1:3, 4);
t3 = A23(1:3, 4);

t11 = simplify(inv(R01)*t1);
t22 = simplify(inv(R12)*t2);
t33 = simplify(inv(R23)*t3);