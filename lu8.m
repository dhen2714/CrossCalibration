function [X]=lu8(AA,BB)
% Solves the problem AX=XB
% using the formulation of
%
% Eight-Space Quaternion Approach 
% for Robotic Hand-eye Calibration.
% Y. Lu, J. C. K. Chou 
%
% Mili Shah
% July 2014

[m,n]=size(AA); n = n/4;

A = zeros(8*n,8);
for i = 1:n
    Ra = AA(1:3,4*i-3:4*i-1);
    Rb = BB(1:3,4*i-3:4*i-1);
    qa = rot2q(Ra); qa = [qa(4); qa(1:3)];
    qb = rot2q(Rb); qb = [qb(4); qb(1:3)];
    ta = [0; AA(1:3,4*i)];
    tb = [0; BB(1:3,4*i)];
    Qc = [qb(1) -qb(2:4)';qb(2:4) qb(1)*eye(3)-skew(qb(2:4))];
    Ql = [qa(1) -qa(2:4)';qa(2:4) qa(1)*eye(3)+skew(qa(2:4))];
    Tc = [tb(1) -tb(2:4)';tb(2:4) tb(1)*eye(3)-skew(tb(2:4))];
    Tl = [ta(1) -ta(2:4)';ta(2:4) ta(1)*eye(3)+skew(ta(2:4))];
    A(8*i-7:8*i-4,:) = [Qc*(Tl-Tc) Ql-Qc];
    A(8*i-3:8*i,1:4) = Ql-Qc;
end
[u,s,v]=svd(A);
v = v(:,8);
nv = norm(v(1:4));
v = v/nv;
R = q2rot(v([2 3 4 1]));
E = [-v(2:4) v(1)*eye(3)+skew(v(2:4))];
t = E*v(5:8);
X = [R t; 0 0 0 1];