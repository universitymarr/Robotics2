function [U, dU_T] = compute_potential_energy_matrix(masses, joint_pos, com_pos, gravity_term, opt_steps)

    % joint_pos: variable of joint position q_i
    % masses: variable of mass of links
    % com_pos: [x, y, z] ( also work with 2D )
    % gravity term: [a, b, g] ( also work with 2D )
    % com_pos and gravity term must have same shape

    n_joint = length(masses);
    U = sym(zeros(1, 1));
    %gravity_term = reshape(gravity_term, length(gravity_term), 1);
    masses = reshape(masses, n_joint, 1);
    joint_pos = reshape(joint_pos, n_joint,1);


    for w = 1:n_joint  
    
        joint_w_mass = masses(w);
        %com_pos_w = com_pos(:,w);
        com_pos_w = com_pos{w};
       
    

        try
            %actual_u = joint_w_mass*transpose(gravity_term)*com_pos_w
            actual_u = joint_w_mass*transpose(gravity_term)*com_pos_w;
            fprintf("Potential energy of link %d: ", w)
            disp(actual_u)

        catch
            disp('Gravity term should be a row vector with shame same of com_i -> need to get a scalare g^T*com_i\n')
        end
        
        U = U + actual_u;
      
    end

    dU_T = transpose(jacobian(U, joint_pos));
    dU_T = collect(dU_T, [joint_pos; gravity_term]);
    U = simplify(U, Steps=opt_steps);
    U = collect(U, [joint_pos; gravity_term]);

    size_dU_t = size(dU_T);

    if size_dU_t(1) > 1 & size_dU_t > 1
        disp("ERROR: vector of COM should contain be contain only z component!")
    end 

end