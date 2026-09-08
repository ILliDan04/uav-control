t0 = 0;     % intial time
t1 = 10;    % end time
eq_num = 2; % number of equations to test
h = 0.001;
hmin = 0.00001;
hmax = 0.1;
emax = 0.01;

eq = {
    struct( ...
        'f', @(t,y) [y(2); -y(1)], ...  % dy1/dt = y2; dy2/dt = -y1;
        'y0', [1; 0], ...
        'true', @(t) [cos(t); -sin(t)] ...
    );

    struct( ...
        'f', @(t,y) [y(1) + y(2); y(2) - y(1)], ...  % dy1/dt = y1 + y2; dy2/dt = -y1 + y2;
        'y0', [1; 0], ...
        'true', @(t) [exp(t).*cos(t); -exp(t).*sin(t)] ...
    );
};

for i = 1:eq_num
    f = eq{i}.f;
    y0 = eq{i}.y0;
    true_res = eq{i}.true;

    res = zeros(size(y0));
    hs = zeros(N,1);
    errors = zeros(N,1);

    [T, Y, H] = rkf45(f, t0, t1, y0, h, emax, hmin, hmax);

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
    legend('RKF45', 'RKF45', 'Exact solution', 'Exact solution');
    grid on;

    hold off;

    % Print relative error value
    err = norm(y1 - true_value(:,end)) / norm(true_value(:,end));
    fprintf("Relative error value for %s: %.15f\n", func2str(f), err);
end