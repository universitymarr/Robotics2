function M = compute_inertia_matrix(kinetic_energy, joint_vel, opt_steps)

    % kinetic_energy: expression of kinetic exergy% kinetic_energy: expression of kinetic exergy
    % joint_vel: variable of joint velocity dq_i
    
    % element is computed as mij = dT/(dqi dqj), double derivative respect
    % to two joints

    T = kinetic_energy;
    n_joint = length(joint_vel);

    M = sym(zeros(n_joint, n_joint));
    joint_vel = reshape(joint_vel, n_joint, 1);

    for i = 1:n_joint
        for j = 1:n_joint
            joint_i = joint_vel(i);
            joint_j = joint_vel(j);
            mij = diff(T, joint_i );
            mij = diff(mij, joint_j);
            M(i,j) = simplify(mij);
        end
    end
    
    M = simplify(M, Steps=opt_steps);

end