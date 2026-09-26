function output = rk4(f, t0, t1, y0, h)
arguments
    f  function_handle % function of t
    t0 (1,1) double    % integration start time
    t1 (1,1) double    % integration end time
    y0 double          % initial solution
    h  (1,1) double    % integration time step
end

if h <= 0
    error("Step size must be positive")
end

if t0 >= t1
    error("Start time must be less than end time")
end

y = y0;
N = ceil((t1 - t0)/h);
t = t0;

for i = 1:N
    h_step = min(h, t1 - t);
    y = rk4_step(f, t, y, h_step);
    t = t + h_step;
end

output = y;

end

function output = rk4_step(f, t, y, h)
arguments
    f function_handle
    t double
    y double
    h double
end
    k1 = f(t,y);
    k2 = f(t + h/2, y + h*k1/2);
    k3 = f(t + h/2, y + h*k2/2);
    k4 = f(t + h, y + h*k3);

    output = y + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

