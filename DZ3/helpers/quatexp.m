function dq = quatexp(dth)
th = norm(dth);
if th < 1e-8
    s = 0.5 - th^2/48;          % sin(th/2)/th, Taylor
    c = 1   - th^2/8;
else
    s = sin(th/2)/th;
    c = cos(th/2);
end
dq = [c; s*dth(:)];