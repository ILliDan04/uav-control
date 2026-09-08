function [T_out, Y_out, H_out] = rkf45(f, t0, t1, y0, h0, emax, hmin, hmax)
arguments
    f       function_handle % function of t
    t0      (1,1) double    % integration start time
    t1      (1,1) double    % integration end time
    y0      double          % initial solution
    h0      (1,1) double    % initial time step
    emax    (1,1) double    % max error value
    hmin    (1,1) double    % min possible value for time step
    hmax    (1,1) double    % max possible value for time step
end

if (emax <= 0)
    error("Error term must be positive")
end

if h0 <= 0 || hmin <= 0 || hmax <= 0
    error("All step size parameters must be positive")
end

if (hmin >= hmax || h0 <= hmin || h0 >= hmax)
    error("Step size bounds are invalid")
end

if t0 >= t1
    error("Start time must be less than end time")
end

T = t0;
Y = y0;
H = [];

t = t0;
y = y0;
h = h0;

while t < t1
    h_next = min(h, t1 - t); % crop the last step
    [y4, y5] = rkf45_step(f, t, y, h_next);

    err = norm(y5 - y4);

    if err > emax && h_next <= hmin
        error("Minimum step size reached before achieving required accuracy");
    end

    if err <= emax
        y = y5;
        Y(:, end + 1) = y;

        t = t + h_next;
        T(end + 1) = t;

        H(end + 1) = h_next;
    end 

    if err == 0
        factor = 5;
    else 
        factor = 0.9 * (emax / err)^(1/5);
        factor = min(5, max(0.2, factor));
    end

    h = factor * h_next;
    h = min(max(h, hmin), hmax);
end

T_out = T;
Y_out = Y;
H_out = H;

end

% ------------------------------------------------------------------------------------------------------

function [y4, y5] = rkf45_step(f, t, y, h)
arguments
    f function_handle
    t double
    y double
    h double
end
k1 = f(t,y);
k2 = f(t + 1/4*h,   y + h*(1/4*k1));
k3 = f(t + 3/8*h,   y + h*(3/32*k1 + 9/32*k2));
k4 = f(t + 12/13*h, y + h*(1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3));
k5 = f(t + h,       y + h*(439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4));
k6 = f(t + 1/2*h,   y + h*(-8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5));

y4 = y + h*(25/216*k1 + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5);
y5 = y + h*(16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6);
end

