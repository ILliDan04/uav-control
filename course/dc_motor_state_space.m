%PARAMETERS

R_a = 2;      % Om
L_a = 0.5;    % N
K_e = 0.1;    % V*s/rad
J_m = 0.01;   % kg * m^2
J_l = 0.02;   % kg * m^2
K_t = 0.1;    % N * m / A
b_   = 0.001; % N * m * s / rad

% parameters symbolic variables
syms Ra Ke La Kt b Jm Jl

% states symbolic variables
syms Ia V w

% EQUATIONS

dIa_dt =    -(Ra/La)    * Ia  -   (Ke/La)   * w + (1/La) * V;
dw_dt  = (Kt/(Jm + Jl)) * Ia  - b/(Jm + Jl) * w;

F = [dIa_dt;dw_dt];
x = [Ia, w];
u = V;

% STATE SPACE FORM

A = jacobian(F, x);
B = jacobian(F, u);

% SUBSTITUTING ACTUAL VALUES

A = double(subs(A, [Ra La Ke Kt b Jm Jl], [R_a L_a K_e K_t b_ J_m J_l]));
B = double(subs(B, [Ra La Ke Kt b Jm Jl], [R_a L_a K_e K_t b_ J_m J_l]));


dc_motor_cfg.ss.A = A;
dc_motor_cfg.ss.B = B;

dc_motor_ss = ss(A, B, [0 1], 0);