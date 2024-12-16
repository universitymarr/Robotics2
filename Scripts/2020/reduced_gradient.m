%% Reduced Gradient
function qd_red = reduced_gradient(JA,JB,T,pd)
    % takes as inputs:
    %   - JA = the jacobian matrix JA of decomposion of J
    %          of size MxM, non singular
    %   - JB = the jacobian matrix JB of decomposion of J
    %          of size Mx(N-M)
    %   - T = identity matrix, corresponding to the q 
    %         variables of the Jacobean JA and JB. 
    %         For example J has dimension 2x3 = [JA JB]. 
    %         JA =[J1 J3] and JB = [J2].
    %         T will be = [1 0 0; 0 0 1; 0 1 0], first column
    %         [1;0;0;] because we have q1, second column is 
    %         [0;0;1;] because we have q3 and the last column is
    %         [0;1;0;] because we have q2
    %   - pd = the cartesian velocity
    %
    % output:
    %   - qd_red = q dot obtain by reduced gradient

    syms q1 q2 q3
    %% JA and JB two matrices give by decomposion that compose J
    J = [JA JB];
   
    %% control if JA is full row rank
    if rank(JA) == size(JA,1)
        invJA=inv(JA);
    else
        disp("my friend, you have a problem")
    end
    
    %% Compute gradient of the objective function given by the problem
    % for example: H(q) = sin^2(q2)+sin^2(q3) 
    % in this case fist row is zero because we don't have variable q1
    gradq_H= [0; 2*sin(q2)*cos(q2); 2*sin(q3)*cos(q3)];
    
    %% Compute reduced gradient
    first_part = invJA*JB;
    transpose_first = -transpose(first_part);
    redgrad_qb_H = [transpose_first 1]*(T*gradq_H);
    
    %% q dot reduced gradient
    qd_red = transpose(T)*[invJA*(pd - JB * redgrad_qb_H); redgrad_qb_H] 
