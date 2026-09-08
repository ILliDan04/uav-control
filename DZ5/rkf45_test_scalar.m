t0 = 0;     % intial time
t1 = 10;    % end time
eq_num = 3; % number of equations to test
h = 0.001;
hmin = 0.00001;
hmax = 0.1;
emax = 0.01;

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

eq_true_res = {
    @(t) exp(2*t);
    @(t) exp(2*t) - t/2 - 1/4;
    @(t) exp(sin(t))
};

for i = 1:eq_num
    f = eq{i};
    y0 = eq_y0(i);
    true_res = eq_true_res{i};

    res = zeros(N,1);
    hs = zeros(N,1);

    [T,Y,H] = rkf45(f, t0, t1, y0, h, emax, hmin, hmax);

    y1 = Y(:, end);
    true_value = true_res(T);

    % Plot y(t)
    figure;
    plot(T, Y, 'o-');
    hold on;
    plot(T, true_value, '-');

    xlabel('t');
    ylabel('y');
    title('RKF45');
    legend('RKF45', 'Exact solution');
    grid on;

    hold off;

    % Print relative error value
    err = abs(y1 - true_value(end)) / abs(true_value(end));
    fprintf("Relative error value for %s: %.15f\n", func2str(f), err);
end