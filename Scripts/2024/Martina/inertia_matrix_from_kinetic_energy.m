function M = inertia_matrix_from_kinetic_energy(T, qdot_vector)
    % Takes as inputs:
    %   - T = the kinetic energy of the robot 
    %   - qdot_vector = the vector of q_dot ex: [qd1,qd2,qd3]
    %
    % Output:
    %   - M = the inertia matrix 
    M = jacobian(jacobian(T, qdot_vector)', qdot_vector);
end