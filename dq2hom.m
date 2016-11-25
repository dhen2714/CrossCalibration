function H=dq2hom(dq)
% Converts a dual quaternion dq into a 4x4 homogeneous matrix H.
% Uses function Qmult.m for quaternion multiplication.

q=dq(:,1);
qe=dq(:,2);
q0=q(1);
qx=q(2);
qy=q(3);
qz=q(4);

R=[(q0^2+qx^2-qy^2-qz^2) 2*(qx*qy-q0*qz) 2*(qx*qz+q0*qy);...
    2*(qx*qy+q0*qz) (q0^2-qx^2+qy^2-qz^2) 2*(qy*qz-q0*qx);...
    2*(qz*qx-q0*qy) 2*(qz*qy+q0*qx) (q0^2-qx^2-qy^2+qz^2)];

q(2:4) = -q(2:4);
t = 2*quatmult(qe,q);
H = [R t(2:4);0 0 0 1];
end
