function u_vector = newton_euler_algorithm(num_joints, DH, qd,qdd, m, I, r, rc, g0, prismatic_indices)
    % This function performs the moving frame algorithm and returns the
    % vectors w, v, vc, and the values of T for each joint
    %
    %inputs:
    % num_joints= an int reapresenting the number of links o fthe robot
    %
    % DH = the DH matrix of the robot in the order [alpha, a, d, theta ]
    %
    % qdts = a vector containing qdots ex: [qd1; qd2; qd3]
    %
    % m
    %
    % I = a vector of kind [I1,I2..]
    %
    % r = a vector containing the values ^ir_{i-1,i} that is the distance
    % between frame i-1 and i computed with respect to frame i 
    %
    % rc = a vector containing the values ^ir_{ci} that is the distance
    % between the i-th CoM from the i-th reference frame
    %
    % g0= a vector representing g wrt fame 0
    %
    % prismatic indices= a list containing the list of prismatic indices,
    % ex for a PPR= [1,2]
    %
    %
    %output:

    syms Ixx Iyy Izz real 


    % Check that the DH is consistent with the indicated number of joints 
    if size(DH, 1) ~= num_joints
        error('Number of rows in DH matrix should be equal to the number of joints.');
    end

    % Check that the r is consistent with the number of joints 
    if size(r,2) ~= num_joints
        size(r,2)
        error('Length of r vector should be equal to the number of joints.');
    end

    % Check that rc is consistent with the number of joints 
    if size(rc,2) ~= num_joints
        error('Length of rc vector should be equal to the number of joints.');
    end

    % Display the number of joints and prismatic indices
    disp(['Number of joints: ', num2str(num_joints)]);

    rotations = sym(zeros(3,3, num_joints));
    rotations_from_zero = sym(zeros(3,3, num_joints));
    i_g=sym(zeros(3,num_joints))

    for i=1:num_joints
        A=DHMatrix(DH(i,:));
        R=A(1:3,1:3);
        rotations(:,:,i)=R;
        A=DHMatrix(DH(1:i,:));
        R=A(1:3,1:3);
        rotations_from_zero(:,:,i)=R;
        i_g(:,i)=rotations_from_zero(:,:,i)'*g0;
    end
    
   

    %forward pass
    % -------- initialization -----------------
    w_prev=[0;0;0];
    w_dot_prev=[0;0;0];
    a_prev=[0;0;0];

    w_vector=zeros(3,num_joints);
    w_dot_vector=zeros(3,num_joints);
    a_vector=zeros(3,num_joints);
    a_c_vector=zeros(3,num_joints);

    % --------- starting forward pass -------------
    for i=1:num_joints
        w=rotations(:,:,i)'*(w_prev+qd(i)*[0;0;1]);
        w=vpa(simplify(w))
        fprintf("the w of link %d is:",i);
        disp(w)
        w_vector(:,i)=w;

        w_dot=rotations(:,:,i)'*(w_dot_prev+qdd(i)*[0;0;1]+cross(qd(i)*w_prev,[0;0;1]));
        fprintf("the w_dot of link %d is:",i);
        disp(w_dot)
        w_dot_vector(:,i)=w_dot;

        a=rotations(:,:,i)'*a_prev+cross(w_dot,r(:,i))+cross(w,r(:,i))
        a=vpa(simplify(a))
        fprintf("the a of link %d is:",i);
        disp(a);
        a_vector(:,i)=a;

        a_c=a+cross(w_dot,r(:,i))+cross(w,cross(w,r(:,i)))
        fprintf("the a_c of link %d is:",i);
        disp(a_c);
        a_c_vector(:,i)=a_c;

        %update the prev
        w_prev=w;
        W_dot_prev=w_dot;
        a_prev=a;
    end

    %starting backward pass
    %---------------- initialization --------------------
    f_next=[0;0;0];
    tau_next=[0;0;0];

    f_vector=sym(zeros(3,num_joints));
    tau_vector=sym(zeros(3,num_joints));

    %--------------- starting backward pass --------------
    for i=num_joints:-1:1
        if i == num_joints
            f=m(i)*(a_c_vector(:,i)-i_g(:,i));
        else
            f=rotations(:,:,i+1)*f_next+m(i)*(a_c_vector(:,i)-i_g(:,i));
        end
        f=vpa(simplify(f));
        fprintf("the f of link %d is:",i);
        disp(f);
        f_vector(:,i)=f;

        if i ==num_joints
            tau=-cross(f,(r(:,i)+rc(:,i)))+I(i)*w_dot_vector(:,i)+cross(w_vector(:,i),(I(i)*w_vector(:,i)));
        else
            tau=rotations(:,:,i+1)*tau_next+ cross((rotations(:,:,i+1)*f_next),rc(:,i))-cross(f,(r(:,i)+rc(:,i)))+I(i)*w_dot_vector(:,i)+cross(w_vector(:,i),(I(i)*w_vector(:,i)));

        end
        tau=vpa(simplify(tau));
        fprintf("the tau of link %d is:",i);
        disp(tau);
        tau_vector(:,i)=tau;
        
        f_next=f;
        tau_next=tau;
    end

    u_vector=sym(zeros(3,num_joints))
    for i=1:num_joints
        if ismember(i, prismatic_indices)
            ui=f_vector(:,i)'*(rotations(:,:,i)'*[0;0;1]);
        else
            ui=tau_vector(:,i)'*(rotations(:,:,i)'*[0;0;1]);
        end
        ui=vpa(simplify(ui));
        fprintf("u of link %d is:",i);
        disp(ui)
        u_vector(:,i)=ui;
    end







