function dq=hom2dq(H)
% Converts a 4x4 homogeneous, rigid body transformation matrix H into a
% dual quaternion represented as a 4x2 matrix dq.
% dq(:,1) corresponds to the real part, and dq(:,2) corresponds to the dual
% part.
% Uses function Qmult.m for quaternion multiplication.

R=[H(1:3,1) H(1:3,2) H(1:3,3)];
t=H(1:3,4);

q0=1+R(1,1)+R(2,2)+R(3,3);
qx=1+R(1,1)-R(2,2)-R(3,3);
qy=1-R(1,1)+R(2,2)-R(3,3);
qz=1-R(1,1)-R(2,2)+R(3,3);
Qcheck=[q0 qx qy qz];

if max(Qcheck)==q0
    q0=sqrt(q0/4);
    qx=(R(3,2)-R(2,3))/(4*q0);
    qy=(R(1,3)-R(3,1))/(4*q0);
    qz=(R(2,1)-R(1,2))/(4*q0);
elseif max(Qcheck)==qx
    qx=-sqrt(qx/4);
    qy=(R(2,1)+R(1,2))/(4*qx);
    qz=(R(1,3)+R(3,1))/(4*qx);
    q0=(R(3,2)-R(2,3))/(4*qx);
elseif max(Qcheck)==qy
    qy=sqrt(qy/4);
    qx=(R(2,1)+R(1,2))/(4*qy);
    qz=(R(3,2)+R(2,3))/(4*qy);
    q0=(R(1,3)-R(3,1))/(4*qy);
elseif max(Qcheck)==qz
    qz=sqrt(qz/4);
    qx=(R(1,3)+R(3,1))/(4*qz);
    qy=(R(3,2)+R(2,3))/(4*qz);
    q0=(R(2,1)-R(1,2))/(4*qz);
end
q=[q0 qx qy qz]';
qprime=0.5*quatmult([0;t],q);
dq=[q qprime];
end