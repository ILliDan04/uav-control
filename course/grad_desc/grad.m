function [KpGrad, KiGrad, KdGrad, TfGrad] = grad(Kp, Ki, Kd, Tf, cost_f, h)
    cur_pid = pid(Kp, Ki, Kd, Tf);
    kp_pid  = pid(Kp + h, Ki, Kd, Tf);
    ki_pid  = pid(Kp, Ki + h, Kd, Tf);
    kd_pid  = pid(Kp, Ki, Kd + h, Tf);
    tf_pid  = pid(Kp, Ki, Kd, Tf + h);

    cur_cost = cost_f(cur_pid);
    cost_kp  = cost_f(kp_pid);
    cost_ki  = cost_f(ki_pid);
    cost_kd  = cost_f(kd_pid);
    cost_tf  = cost_f(tf_pid);

    KpGrad = (cost_kp - cur_cost) / h;
    KiGrad = (cost_ki - cur_cost) / h;
    KdGrad = (cost_kd - cur_cost) / h;
    TfGrad = (cost_tf - cur_cost) / h;
end