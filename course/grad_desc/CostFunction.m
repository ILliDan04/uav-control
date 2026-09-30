function cost = CostFunction(Q, R, wRT, sys_plant, sys_pid, T, h, r)

    t = 0:h:T;
    Ref = r * ones(size(t));

    t = t(:);
    Ref = Ref(:);

    late = t >= 2;

    T = feedback(sys_plant * sys_pid, 1);

    yOut = lsim(T, Ref, t);
    yOut = yOut(:);

    recover_u = feedback(sys_pid, sys_plant);

    uIn = lsim(recover_u, Ref, t);
    uIn = uIn(:);

    e = Ref - yOut;  

    J_error = trapz(t, sqrt(e.^2 + 0.01^2) - 0.01);
    J_effort = trapz(t, uIn.^2);
    J_rise_time = mean(max(0, 1 - yOut(late)/r).^2);

    cost = Q * J_error + R * J_effort + wRT * J_rise_time;
end