function q=quatmult(q1,q2)
% Multiplies two quaternions q1 and q2 represented by 4x1 vectors.
% Output q is a 4x1 vector.

q1vec=q1(2:4);
q2vec=q2(2:4);
q1s=q1(1);
q2s=q2(1);

qs=q1s*q2s-q1vec'*q2vec;
qvec=(q1s*q2vec+q2s*q1vec+cross(q1vec,q2vec))';

q=[qs qvec]';

