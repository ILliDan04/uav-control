function c = quatmul(a, b)      % Hamilton product, scalar first
aw = a(1); av = a(2:4);
bw = b(1); bv = b(2:4);
c = [aw*bw - av'*bv;
     aw*bv + bw*av + cross(av, bv)];