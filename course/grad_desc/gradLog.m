function g = gradLog(th, cost_f, d)
% th = [log Kp; log Ki; log Kd; log Tf]
    g = zeros(numel(th),1);
    for i = 1:numel(th)
        e = zeros(numel(th),1); e(i) = d;
        g(i) = (cost_f(mkpid(th+e)) - cost_f(mkpid(th-e))) / (2*d);
    end
end