t0 = 0;     % intial time
t1 = 10;    % end time
eq_num = 3; % number of equations to test
dt = 0.0001;
N = 1000;   % max number of steps

eq = {
    % dy/dt=2y -> y=e^2t
    @(t, y) 2*y;

    % dy/dt=2y -> y=e^2t - t/2 - 1/4
    @(t, y) 2*y + t;

    % dy/dt=ycos(t) -> y=e^sin(t)
    @(t, y) y*cos(t)
};

eq_y0 = [
    1;
    3/4;
    1
];

eq_true_res = [
    exp(2*t1);
    exp(2*t1) - t1/2 - 1/4;
    exp(sin(t1))
];

for i = 1:eq_num
    f = eq{i};
    y0 = eq_y0(i);
    true_res = eq_true_res(i);

    res = zeros(N,1);
    hs = zeros(N,1);

    for n = 1:N 
        h = n * dt;
        hs(n) = h;

        res(n) = rk4(f, t0, t1, y0, h);
    end

    errors = abs(res - true_res) / true_res;

    % Plot error
    figure;
    plot(hs, errors, 'o-');

    xlabel('Step size h');
    ylabel('Relative error');
    title('RK4 error vs. step size');
    grid on;
end