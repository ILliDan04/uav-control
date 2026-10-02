function cost = CostFunction(Q, R, wRT, wPk, sys_plant, sys_pid, Tend, h, r, uMax)
    t   = (0:h:Tend)';
    Ref = r*ones(size(t));
    late = t >= 2;

    Tcl = feedback(sys_plant*sys_pid, 1);
    Gu  = feedback(sys_pid, sys_plant);        % u/r = C/(1+CP)

    if any(real(pole(Tcl)) > 0)                
        cost = 1e3; info = struct('unstable',true); return
    end

    yOut = lsim(Tcl, Ref, t);  yOut = yOut(:);
    uIn  = lsim(Gu,  Ref, t);  uIn  = uIn(:);

    eN = (Ref - yOut)/r;
    uN = uIn/uMax;
    d  = 0.01;

    J_error  = mean(sqrt(eN.^2 + d^2) - d);
    J_effort = mean(uN.^2);
    J_rise   = mean(max(0, 0.9 - yOut(late)/r).^2);
    J_peak   = mean(max(0, abs(uN) - 0.95).^2);

    cost = Q*J_error + R*J_effort + wRT*J_rise + wPk*J_peak;
end