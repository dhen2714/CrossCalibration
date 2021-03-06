function X=dqcrosscalib(A,B)
% X is the cross calibration matrix AX=XB between two coordinate systems.
% Ai,Bi are 4x4 matrices representing rigid body motion in the two frames.
% The inputs of this function A and B are 4x4xN matrices representing a 
% series of N motions in the two frames.

n=size(A,3);
T=zeros(8,6*n);

for i=1:n
    Ai=A(:,:,i);
    Bi=B(:,:,i);
    a=hom2dq(Ai);
    b=hom2dq(Bi);
    S=[(a(2:4,1)-b(2:4,1)),skew(a(2:4,1)+b(2:4,1)),zeros(3,4);...
       (a(2:4,2)-b(2:4,2)),skew(a(2:4,2)+b(2:4,2)),(a(2:4,1)-b(2:4,1)),skew(a(2:4,1)+b(2:4,1))];
    T(:,((i-1)*6+1):(i*6))=S';
end
[u, s, v]=svd(T');

u1=v(1:4,7);
v1=v(5:8,7);
u2=v(1:4,8);
v2=v(5:8,8);

a=u1'*v1;
b=u1'*v2+u2'*v1;
c=u2'*v2;

s=roots([a b c]);
val1=s(1)^2*(u1'*u1)+2*s(1)*(u1'*u2)+u2'*u2;
val2=s(2)^2*(u1'*u1)+2*s(2)*(u1'*u2)+u2'*u2;

if val1>val2
    s=s(1);
    val=val1;
else
    s=s(2);
    val=val2;
end

l2=sqrt(1/val);
l1=s*l2;
q=l1*v(:,7)+l2*v(:,8);
X=dq2hom([q(1:4) q(5:8)]);
end

    