function [C, c_sub] = compute_christoffel_matrix(inertia_matrix, joint_pos, joint_vel, opt_steps)
    
    % inertia_matrix: expression of inertia matrix
    % joint_pos: variable of joint position q_i
    % ch_matrix_disp: can be true or false, it is used to display the
    % single christoffel matrix

    c_sub = {};

    M = inertia_matrix;
    n_joint = length(joint_pos);
    C = sym(zeros(n_joint, 1));
    joint_pos = reshape(joint_pos, n_joint, 1);
    joint_vel = reshape(joint_vel, n_joint, 1);

    for k = 1:n_joint
        joint = joint_pos(k);
        mk = M(:,k); % inertia matrix of joint k, its a columns
    
        %term1 = jacobian(mk,joint_pos);
        %term2 = jacobian(mk,joint_pos);
        %term3 = diff(M,joint);
        %ck = (1/2)*(term1 + transpose(term2) - term3 );

        ck = (1/2)*(jacobian(mk,joint_pos)+transpose(jacobian(mk,joint_pos))-diff(M,joint));
        ck = simplify(ck, Steps=opt_steps);
        c_sub{k} = ck;

        C(k) = simplify(transpose(joint_vel)*ck*joint_vel);
    end
    
    C = simplify(C, Steps=opt_steps);
    C = collect(C, joint_vel);
end