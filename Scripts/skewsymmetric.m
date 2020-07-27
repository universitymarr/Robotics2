
%% Skewsymmetric
syms a1 a2 a3 a4 real
syms q1 q2 q3 real
syms dq1 dq2 dq3 real
%% inertia matrix
q = [q1;q2]

M = [ a1+a2*sin(q2)^2+a3*cos(q2)^2 0;
      0 a4]

%% coriolis
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1))
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2))

pause

dq=[dq1;dq2];
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c=simplify([c1;c2])

pause
%% M dot

dM = diff(M,q1)*dq1+diff(M,q2)*dq2

pause
%% skewsymm

S1 = [dq'*C1;
      dq'*C2;]
 
pause
disp("skewsymmetric matrix")

skew= simplify(dM-2*S1)
pause
%% not skewsymm

S2= [dq2*sin(2*q2)*(a2-a3) 0;
     -dq1*sin(2*q2)*(a2-a3) 0]

disp(" not skewsymmetric matrix")
not_skew=simplify(dM-2*S2)
