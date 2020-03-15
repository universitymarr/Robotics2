%% Jacobian Based Method - Weighted Pseudoinverse
function [qwd, JWpse] = JBM_weightedpseudo(J,pd)
    % takes as inputs:
    %   - J = the jacobian matrix
    %   - pd = the cartesian velocity
    %
    % output:
    %   - qwd = weighted q dot
    %   - JWpse = weighted pseudoinverse 


    %% define W, weighted matrix
    syms w
    W=[1 0 0; 0 1 0; 0 0 w];

    %% Control if J is full row rank
    if rank(J) == size(J,1)
        % compute Weighted Pseudoinverse
        secondPart = simplify(inv(J*inv(W)*transpose(J)));
        JWpse = simplify(inv(W)*transpose(J)*secondPart);
        pause
        detJW = simplify(det(J*inv(W)*transpose(J)));
        pause
    else
        disp('Oh my friend, there is a problem')
    end

    %% general solution to minimize objective function and understand 
    %% which value we can choose for w
    vecdot= [sym('qd1'); sym('qd2'); sym('qd3');]; 

    H = (1/2)*transpose(vecdot)*W*vecdot
    pause

    %% compute joint velocity
    qwd = simplify(JWpse * pd);

