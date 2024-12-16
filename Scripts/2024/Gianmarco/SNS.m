function [new_values, old_values]= SNS(task, J_array, constraint_min, constraint_max, velocity_initial_condition)
    % This algorithm works for velocity and acceleration.
    % Change only equation to calculate desired task

    % PARAMETERS:
    %   - Task: task to achieve
    %   - J_array: array which contain the jacobian matrix. It is:
    %       - {J}       : For velocity
    %       - {J, dJ}   : For acceleration
    %   - constrain_min/max: Array containing constraints for each joint
    %   - velocity_initial_condition: Needed only for acceleration case, else you can put to 0 -> DEPRECATED

    % RETURNS:
    %   - new_values: Array containing the new values for the joints
    %   - old_values: Array containing the old values for the joints

    n_joint = length(constraint_min);
    task = reshape(task, length(task), 1);
    constraint_min = reshape(constraint_min, n_joint, 1);
    constraint_max = reshape(constraint_max, n_joint, 1);

    init_value = -999;
    new_values = linspace(init_value, init_value, n_joint);

    % Reshape the velocity_initial_condition to a column vector [num_joints x 1]
    if length(J_array) == 2
        velocity_initial_condition = reshape(velocity_initial_condition, n_joint, 1);
    end

    % Check if J_array is {J} or {J, dJ}
    % This means that we are working with velocity or acceleration respectively
    if isscalar(J_array)
        disp("---- SNS for velocity")
        J = J_array{1};
        type_task = 'v';
    elseif length(J_array) == 2
        disp("---- SNS for acceleration")
        J = J_array{1};
        dJ = J_array{2};
        type_task = 'a';
    else
        disp("Error: J_array = {J, dJ} or J_array = {J}")
        return
    end

    % Check based on the type of task
    if type_task == 'v'
        old_values = pinv(J)*task;
    elseif type_task == 'a'
        old_values = pinv(J)*(task - dJ*velocity_initial_condition);
    end

    for i = 1:n_joint
        fprintf("==========> Loop number: %d", i)
        if type_task == 'v'
            result = pinv(J)*task;
        elseif type_task == 'a'
            result = pinv(J)*(task - dJ*velocity_initial_condition);
        end

        result = reshape(result, n_joint-i+1, 1);
        disp("Result: ")
        disp(result)

        % picking maximum violating constraint for max
        violating_max_logical = logical(result > constraint_max);
        violating_max_value = max(dot(result - constraint_max, violating_max_logical));
        
        % picking maximum violating constraint for min
        violating_min_logical = logical(result < constraint_min);
        violating_min_value = min(dot(result - constraint_min, violating_min_logical));

        % Print the violating constraints
        fprintf("Violating min: %d\n", violating_min_value)
        fprintf("Violating max: %d\n", violating_max_value)

        if all(~violating_min_logical) && all(~violating_max_logical)
            logic_subs = logical(new_values == init_value);
            
            counter_subs = 1;
            % Substitute the value in the new_values array
            for j = 1:length(new_values)
                if logic_subs(j) == 1
                    new_values(j) = result(counter_subs);
                    counter_subs = counter_subs + 1;
                end
                
            end
            
            % Reshape the new_values array to a column vector [num_joints x 1]
            new_values = reshape(new_values, n_joint, 1);
            break
        end
        
        % Check which constraint is violating the most
        if abs(violating_min_value) > abs(violating_max_value)
            violating_joint_index = find( (result-constraint_min) == violating_min_value );
            violation_type = 'violating_min';
        else
            violating_joint_index = find( (result-constraint_max) == violating_max_value );
            violation_type = 'violating_max';
        end

        fprintf("--- Joint %d exceeded bounds ==> %f\n", violating_joint_index, result(violating_joint_index));
        
        if strcmp(violation_type, 'violating_min')
            new_values(violating_joint_index) = constraint_min(violating_joint_index);
        else
            new_values(violating_joint_index) = constraint_max(violating_joint_index);
        end
        fprintf("--- Joint %d now set at %f \n", violating_joint_index, new_values(violating_joint_index))

        if type_task == 'v'
            selected_J = J(:, violating_joint_index);
            actual_velocity = new_values(violating_joint_index);
            task = task - selected_J*actual_velocity;

        elseif type_task == 'a'
            selected_J = J(:, violating_joint_index);

            actual_acceleration = new_values(violating_joint_index);

            task = task - selected_J*actual_acceleration;
 
            dJ(:,violating_joint_index) = [];
            velocity_initial_condition(violating_joint_index) = [];
        end
        
        % Empty all the values related to the violating joint
        J(:, violating_joint_index) = [];
        constraint_max(violating_joint_index) = [];
        constraint_min(violating_joint_index) = [];
     end
end