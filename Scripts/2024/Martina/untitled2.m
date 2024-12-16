clear all
clc
syms q1 q2 q3 l y real
digits(4);
%%
l=1

p=[l*(cos(q1)+cos(q1+q2)+cos(q1+q2+q3));
    l*(sin(q1)+sin(q1+q2)+sin(q1+q2+q3))]
j=jacobian(p,[q1,q2,q3])
jp=pinv(j)
%vpasolve((y-2)^2+(0.25^2-0.5^2),y)
H=norm([cos(q1)+cos(q1+q2)-1;sin(q1)+sin(q1+q2)-1.567])
DeH=vpa(subs(simplify([diff(H,q1);diff(H,q2);diff(H,q3)]),[q1,q2,q3],[0,pi/2,-pi/2]))
first_part=jp*[0;1]
second_part=eye(3)-jp*j
qdot=vpa(subs(first_part+second_part*-DeH,[q1,q2,q3],[0,pi/2,-pi/2]))
vm=subs(j(:,1:2)*qdot(1:2),[q1,q2,q3],[0,pi/2,-pi/2])
qdot1=vpa(subs(jp*[0;1],[q1,q2,q3],[0,pi/2,-pi/2]))
P1=vpa(simplify(subs(-jp*j,[q1,q2,q3],[0,pi/2,-pi/2])))
p2=[l*(cos(q1)+cos(q1+q2));
    l*(sin(q1)+sin(q1+q2))]
j2=subs(jacobian(p,[q1,q2,q3]),[q1,q2,q3],[0,pi/2,-pi/2])
qdot2=qdot1+pinv((j2*P1))*(vm-j2*qdot1)