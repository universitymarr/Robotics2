num_joints = 3; % Number of joints
m = 2; % Task dimension
V = [-1 , -1.5 , -2;
      1 ,  1.5 , 2]; % First row lower limits, second row upper limits
J = [-1   , -1    , -0.5  ;
    -0.366, -0.866, -0.866]; % Jacobian

task = [2; 1]; % Task

% SNS (task, J_array, constraint_min, constraint_max, velocity_initial_condition)
[new_values, old_values] = SNS(task, {J}, V(1,:)', V(2,:)', zeros(num_joints, 1));