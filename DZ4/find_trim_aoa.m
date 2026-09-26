f = @(x) 2.7939 * 9.80665 * cos(x) - 8.2 * 60.8837 * x;
x0 = 0.01;
x = fzero(f, x0);
disp(x);
disp(f(x));
disp(90 * cos(x));
disp(90 * sin(x));