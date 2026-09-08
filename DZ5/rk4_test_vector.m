t0 = 0;     % intial time
t1 = 10;    % end time
eq_num = 2; % number of equations to test
dt = 0.0001;
N = 1000;   % max number of steps

eq = {
    struct( ...
        'f', @(t,y) [y(2); -y(1)], ...  % dy1/dt = y2; dy2/dt = -y1;
        'y0', [1; 0], ...
        'true', @(t) [cos(t); -sin(t)] ...
    );

    struct( ...
        'f', @(t,y) [y(1) + y(2); y(2) - y(1)], ...  % dy1/dt = y1 + y2; dy2/dt = -y1 + y2;
        'y0', [1; 0], ...
        'true', @(t) [exp(t)*cos(t); -exp(t)*sin(t)] ...
    );
};

for i = 1:eq_num
    f = eq{i}.f;
    y0 = eq{i}.y0;
    true_res = eq{i}.true(t1);

    res = zeros(size(y0));
    hs = zeros(N,1);
    errors = zeros(N,1);

    for n = 1:N 
        h = n * dt;
        hs(n) = h;

        res = rk4(f, t0, t1, y0, h);

        errors(n) = norm(res - true_res) / norm(true_res);
    end

    % Plot error
    figure;
    plot(hs, errors, 'o-');

    xlabel('Step size h');
    ylabel('Relative error');
    title('RK4 error vs. step size');
    grid on;
end