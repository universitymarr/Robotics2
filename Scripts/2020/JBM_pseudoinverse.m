%% Jacobian Based Method - Pseudoinverse
function [qd, Jpse] = JBM_pseudoinverse(J,pd)
    % takes as inputs:
    %   - J = the jacobian matrix
    %   - pd = the cartesian velocity
    %
    % output:
    %   - qd = q dot 
    %   - Jpse = pseudoinverse 

    %% Control if J is full row rank
    if rank(J) == size(J,1)
        % compute Pseudoinverse
        Jpse = simplify(transpose(J)*(inv(J*(transpose(J)))))
        detJ= simplify(det(J*(transpose(J))))
    else
        Jpse = pinv(J)     
    end

    %% compute joint velocity
    qd = simplify(Jpse * pd)
