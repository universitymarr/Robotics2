function [w, v, vc, T] = moving_frames_algorithm(num_joints, DH, qdots, m,  rc, prismatic_indices,I)
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
    % m = a vector of masses [m1 ; ... ; mn]
    %
    %
    % rc = a vector containing the values ^ir_{ci} that is the distance
    % between the i-th CoM from the i-th reference frame
    %
    % prismatic indices= a list containing the list of prismatic indices,
    % ex for a PPR= [1,2]
    %
    % I = a vector of 3 by 3 matrices [I1,..,In]
    %
    %output:

    syms Ixx Iyy Izz real
    


    % Check that the DH is consistent with the indicated number of joints 
    if size(DH, 1) ~= num_joints
        error('Number of rows in DH matrix should be equal to the number of joints.');
    end


    % Check that rc is consistent with the number of joints 
    if size(rc,2) ~= num_joints
        error('Length of rc vector should be equal to the number of joints.');
    end

    % Display the number of joints and prismatic indices
    disp(['Number of joints: ', num2str(num_joints)]);
    

    %START THE ALGORITHM:
    v00=[0;0;0]
    w00=[0;0;0]

    % Initialize outputs
    w = zeros(3, num_joints);
    v = zeros(3, num_joints);
    vc = zeros(3, num_joints);
    T = cell(num_joints);
    
    
    %Start the loop
    for i = 1:num_joints
        In=I(1:3,(i-1)*3+1:i*3)
      
        
        if ismember(i, prismatic_indices)
            sigma = 1;
        else
            sigma =0;
        end
        A=DHMatrix(DH(i,:));
        R=A(1:3,1:3)
        p=A(1:3,4)

        r_i=R'*p
        % computing omega
        if i ==1
            w_prev= w00;
        else
            w_prev= wi
        end
        
        wi=R'*(w_prev + (1- sigma)*[0;0;qdots(i)]);
        fprintf('the value of w_%d is:',i);
        disp(wi);
        
        % computing v
        if i ==1
            v_prev= v00;
        else
            v_prev= vi;
        end 
        
        vi=R'*(v_prev + sigma*[0;0;qdots(i)])+ cross(wi,r_i,1);
        fprintf('The value of v_%d is:', i);
        disp(vi)

        %computing vc
        vci=vi+cross(wi,rc(:,i));
        fprintf('The value of vc_%d is:',i);
        disp(vci)

        %computing T
        Ti=simplify(0.5*m(i)*(vci'*vci)+0.5*wi'*In*wi);
        fprintf('The value of T_%d is:', i);
        disp(Ti)
        T{i}=Ti;
    
    
end