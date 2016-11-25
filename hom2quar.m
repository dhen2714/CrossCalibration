function dq = hom2quar(H)
% Converts 4x4 homogeneous matrix H
% to a dual quaternion dq represented as
% dq(:,1) + e*dq(:,2)
%
% Mili Shah

R = H(1:3,1:3);
t = H(1:3,4);

R = logm(R);
r = [R(3,2) R(1,3) R(2,1)]';
theta = norm(r);
l = r/norm(theta);

q = [cos(theta/2); sin(theta/2)*l];
qprime = .5*Qmult([0;t],q);
dq=[q qprime];
end