function [c, C] = inertia_matrix_to_coriolis(M, q, dq)
    % Takes as inputs:
    %   - M = the inertia matrix
    %   - q = a vertical vector of q values
    %   - dq = a vertical vector dot_q
    % Output:
    %   - c = robot centrifugal and Coriolis term
    %   - C = Christoffel matrices

    % Initialize variables
    n = length(q); % Number of joints
    C = cell(n, 1); % Initialize cell array to store Christoffel matrices
    
    % Loop over each joint
    for i = 1:n
        % Compute Christoffel matrix for the i-th joint
        disp(['Christoffel matrix for joint ', num2str(i)])
        Mi = M(:, i); % Select the i-th column of the inertia matrix
        Ci = (1/2) * (jacobian(Mi, q) + jacobian(Mi, q)' - diff(M, q(i)));
        C{i} = Ci; % Store the Christoffel matrix
        
        % Display the Christoffel matrix
        disp(['C', num2str(i), ' = ']);
        disp(Ci);
    end
    
    % Compute robot centrifugal and Coriolis terms
    disp("Robot centrifugal and Coriolis terms")
    c = sym(zeros(n, 1));
    for i = 1:n
        c(i) = dq' * C{i} * dq;
        disp(['c', num2str(i), ' = ']);
        disp(c(i));
    end
end
