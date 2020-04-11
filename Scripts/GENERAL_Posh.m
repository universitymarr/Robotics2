%% PRP-planar robot (Midterm 2018-2019)

% 1) N.B: --> each column refers to a sigle joint's parameters
% 2) GAME STARTS

% number of joints 
N = 3;

% 'q' and 'dq' are respectively vectors of joint's variables and their derivatives
q = sym(['q1';'q2';'q3']);           assume(q, 'real');           
dq = sym(['dq1';'dq2';'dq3']);       assume(dq, 'real');  

% gravity vector
g0 = [0; 0; 0];             % Horizontal plane otherwise g0=[0; -sym('g0'); 0];

% constants
k1=sym('k1');              assume(k1, 'real');   
D2=sym('D2');              assume(D2, 'real');
L2=sym('L2');              assume(L2, 'real');

% masses
m = sym(['m1'; 'm2'; 'm3']);     assume(m, 'real');        


% Inertia 
Iz = [sym('Iz1'); sym('Iz2'); sym('Iz3')];  assume(Iz, 'real');
I = cell(1,N);

for i=1:N
    I{i} = [0  0   0;                     
            0  0   0;
            0  0   Iz(i)];
end
     
% positions of CoM's w.r.t. the base frame  
rc =  [0      D2*cos(q(2))            L2*cos(q(2))-q(3)*sin(q(2));           
       q(1)   k1+q(1)+D2*sin(q(2))    q(1)+k1+L2*sin(q(2))+q(3)*cos(q(2));
       0             0                            0];       
   
% angular velocities
w = [0   0       0;
     0   0       0;
     0  dq(2)  dq(2)];    
 
% traslational velocities 
Vc = [  0          -D2*dq(2)*sin(q(2))                -L2*dq(2)*sin(q(2))-(dq(3)*sin(q(2))+q(3)*dq(2)*cos(q(2)));
        dq(1)      dq(1)+D2*dq(2)*cos(q(2))         dq(1)+L2*dq(2)*cos(q(2))+(dq(3)*cos(q(2))-q(3)*dq(2)*sin(q(2)));
        0                     0                                           0];                                            
  
%% Kinetic Energy and Inertia Matrix    
T = sym(zeros(N,1));
Tot_Kin_En = 0;

for i=1:N
    T(i) = 0.5*m(i)*(Vc(:,i)'*Vc(:,i))+ 0.5*w(:,i)'*I{i}*w(:,i);
    T(i) = simplify(T(i),'IgnoreAnalyticConstraints',true);
    Tot_Kin_En = Tot_Kin_En + T(i);
end

Tot_Kin_En = simplify(collect(Tot_Kin_En,dq));

M = sym(zeros(N,N)); 
for i=1:N
    for j=1:N
        if i==j
            Tot_Kin_En = collect(Tot_Kin_En, dq(i)^2);
            M(i,i) = simplify(diff(Tot_Kin_En, dq(i), 2));
       
        else
            Tot_Kin_En = collect(Tot_Kin_En, dq(i)*dq(j));
            tmp = simplify(diff(Tot_Kin_En, dq(i)));
            M(i,j) = simplify(diff(tmp, dq(j)));
            M(j,i) = M(i,j);
        end             
    end
end
M=simplify(M);

%% Calculate Christoffel matrix
c = sym(zeros(N,1));
C = cell(1,N);

for j=1:N
    C{j} = 0.5*simplify((jacobian(M(:,j),q) + transpose(jacobian(M(:,j),q)) - diff(M,q(j)))); 
    c(j) = dq'*C{j}*dq;
    c(j) = simplify(collect(c(j), dq));
end

%% Calculate Gravitational term
U = 0;
aux = sym(zeros(N,1));                                                                
                                                                                                  
for i=1:N   
    aux(i) = simplify(-m(i)*g0'*rc(:,i)); 
    U = U + aux(i);
end

G = jacobian(U,q)';

% GAME OVER 


















