%% Jacobian Based Method - Weighted Pseudoinverse

%% IMPORTANT : BEFORE RUNNING THIS CODE, RUN THE FILE OF THE CORRESPONDIGN ROBOT
%% FOR EXAMPLE : planar_3R for 3R planar robot

%% define W, weighted matrix

% if we have joints with diffent nature, we can consider W in such way
% for example : PPR robot, we have 2 joints that are prismatic and 1
% revolut so we can consider 2 first joints with same weight and revolut
% proportional to other two

syms w
W=[1 0 0; 0 1 0; 0 0 w]

%% Control if J is full row rank
if rho == size(J,1)
    % compute Weighted Pseudoinverse
    secondPart = simplify(inv(J*inv(W)*transpose(J)));
    JWpse = simplify(inv(W)*transpose(J)*secondPart)
    pause
    detJW = simplify(det(J*inv(W)*transpose(J)))
    pause
else
    disp('Oh my friend, there is a problem')
end

%% general solution to minimize objective function and understand 
%% which value we can choose for w
vecdot= [sym('qd1'); sym('qd2'); sym('qd3');] 

H = (1/2)*transpose(vecdot)*W*vecdot
pause
% to choose in this case can observe determinant detJW = 1+l^2/w
% we choose w = l^2 to obtain determinant equal to 2

W =[1 0 0; 0 1 0; 0 0 l^2];
secondPart = inv(J*inv(W)*transpose(J));
JWpse = simplify(inv(W)*transpose(J)*secondPart)

%% compute joint velocity
qwd = simplify(JWpse * pd)
pause
% joint velocity with values
qwd_val= simplify(subs(qwd,{l,q3},{0.5,pi/6}))

%% NB 
% if we choose a value of w that is too unbalanced compared to the others 
% we end up penalizing the movement of that joint.
