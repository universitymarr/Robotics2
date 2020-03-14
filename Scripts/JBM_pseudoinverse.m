%% Jacobian Based Method - Pseudoinverse

%% IMPORTANT : BEFORE RUNNING THIS CODE, RUN THE FILE OF THE CORRESPONDIGN ROBOT
%% FOR EXAMPLE : planar_3R for 3R planar robot

%% Control if J is full row rank
if rho == size(J,1)
    % compute Pseudoinverse
    Jpse = simplify(transpose(J)*(inv(J*(transpose(J)))))
    detJ= simplify(det(J*(transpose(J))))
else
    Jpse = pinv(J)     
end

%% compute joint velocity
qd = simplify(Jpse * pd)
pause

% joint velocity with values
qd_val= simplify(subs(qd,{l,q3},{0.5,pi/6}))
