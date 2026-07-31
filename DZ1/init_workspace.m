clear P;

%% ------------------------------------------------------------------
%  Environment
%  ------------------------------------------------------------------
P.g = 9.81;                 % [m/s^2]  acceleration due to gravity

%% ------------------------------------------------------------------
%  Mass properties  (thesis section 4.2, eq. 4.7, Table 4.2)
%  ------------------------------------------------------------------
% Moment of inertia tensor of the whole quadrotor, diagonal because the
% body axes coincide with the principal axes (eq. 2.21 / 4.7)
P.Ix = 0.0093;              % [kg*m^2] roll  inertia (about body x)
P.Iy = 0.0092;              % [kg*m^2] pitch inertia (about body y)
P.Iz = 0.0151;              % [kg*m^2] yaw   inertia (about body z)
P.I  = diag([P.Ix P.Iy P.Iz]);

% Total mass. The thesis never states it as a single number; it is
% recovered from the hover condition used to identify b (eq. 4.9):
%   b = m*g / (4*Omega0^2)  =>  m = 4*b*Omega0^2 / g
% With b = 1.5108e-5 and Omega0 = 76 RPS this gives ~1.40 kg, consistent
% with the sum of the part masses in Table 4.2 (~1.31 kg + props/wiring).
P.Omega0 = 463.1;       % [rad/s]  rotor speed at hover trim (76 RPS)
P.b      = 1.5108e-5;       % [N*s^2/rad^2 = kg*m] propeller thrust coeff.
                            %   T_i = b * Omega_i^2          (eq. 4.8/4.9)
P.m      = 4*P.b*P.Omega0^2 / P.g;   % [kg] total mass 

%% ------------------------------------------------------------------
%  Geometry
%  ------------------------------------------------------------------
P.l = 0.214;                % [m] lever arm: distance motor axis -> CG
                            %   (Table 4.2: motors at 21.4 cm)

%% ------------------------------------------------------------------
%  Propeller aerodynamics  (thesis section 4.3)
%  ------------------------------------------------------------------
% Thrust:      T_i  = b * Omega_i^2                (P.b defined above)
% Drag torque: Q_i  = d * Omega_i^2                (eq. 4.19)
P.d = 4.406e-7;             % [N*m*s^2/rad^2 = kg*m^2] propeller drag factor

%% ------------------------------------------------------------------
%  Rotating masses (gyroscopic terms, thesis sections 2.3.1 and 4.1-4.3)
%  ------------------------------------------------------------------
P.Jm = 2.506e-6;            % [kg*m^2] motor rotor inertia (eq. 4.22)
P.Ip = 4.439e-5;            % [kg*m^2] propeller inertia   (below Table 4.1)
P.Jr   = P.Jm + P.Ip;         % [kg*m^2] total spinning inertia per rotor
                            %   (eq. 4.14: J_r = J_m + I_p)

% Total z-axis angular momentum of the four rotors (eqs. 2.29-2.31):
%   H_z = Jr * (s1*Omega1 + s2*Omega2 + s3*Omega3 + s4*Omega4)
% where s_i = +/-1 encodes each rotor's spin direction about body z.
% Pairs (1,3) and (2,4) counter-rotate
P.rotor_spin_sign = [ -1 1  -1 1 ];   % [-] s_i for rotors 1..4

%% ------------------------------------------------------------------
%  Motor (engine) dynamics  (thesis sections 2.3.2, 4.1.4, eq. 4.6)
%  ------------------------------------------------------------------
% Identified first-order model, PPM input -> rotor speed [RPS], valid
% around the hover trim point:
%   G_m(s) = K / (tau*s + 1)     <=>   wdot_m = -(1/tau)*w_m + (K/tau)*u
P.motor.K   = 0.7;          % [-] DC gain (PPM units -> RPS)
P.motor.tau = 0.1;          % [s] time constant
P.motor.tau_np = 0.5;       % [s] time constant WITHOUT propeller (eq. 4.17
                            %     measurement) - for reference only

% Rotor speed saturation used in the thesis Simulink model (chapter 5)
P.motor.rps_min = 0;        % [RPS] lower limit
P.motor.rps_max = 150;      % [RPS] upper limit
P.motor.omega_min = P.motor.rps_min * 2*pi;   % [rad/s]
P.motor.omega_max = P.motor.rps_max * 2*pi;   % [rad/s]

%% ------------------------------------------------------------------
%  Derived quantities useful for control design
%  ------------------------------------------------------------------
P.hover.T_total  = P.m * P.g;             % [N]  total hover thrust
P.hover.T_rotor  = P.hover.T_total / 4;   % [N]  per-rotor hover thrust
P.hover.Hz       = P.Jr * (P.rotor_spin_sign * ...
                   (P.Omega0*ones(4,1)));  % [kg*m^2/s] residual rotor
                                           % momentum at hover (~0 if the
                                           % speeds are matched)

% Control effectiveness at hover (useful for linearization):
% dT/dOmega = 2*b*Omega0, dMx/dOmega = 2*l*b*Omega0, dMz/dOmega = 2*d*Omega0
P.hover.dT_dOmega  = 2*P.b*P.Omega0;      % [N*s/rad]   per rotor
P.hover.dMx_dOmega = 2*P.l*P.b*P.Omega0;  % [N*m*s/rad] per roll-pair rotor
P.hover.dMz_dOmega = 2*P.d*P.Omega0;      % [N*m*s/rad] per rotor

%% ------------------------------------------------------------------
%  Export to base workspace as individual variables (optional)
%  ------------------------------------------------------------------
g  = P.g;    m  = P.m;
Ix = P.Ix;   Iy = P.Iy;   Iz = P.Iz;   I = P.I;
l  = P.l;    b  = P.b;    d  = P.d;
Jm = P.Jm;   Ip = P.Ip;   Jr = P.Jr;
Omega0 = P.Omega0;

fprintf('Quadrotor parameters loaded:\n');
fprintf('  m  = %.4f kg,  l = %.3f m\n', P.m, P.l);
fprintf('  I  = diag([%.4f %.4f %.4f]) kg*m^2\n', P.Ix, P.Iy, P.Iz);
fprintf('  b  = %.4e kg*m,  d = %.4e kg*m^2\n', P.b, P.d);
fprintf('  Jr = %.4e kg*m^2 (Jm = %.3e + Ip = %.3e)\n', P.Jr, P.Jm, P.Ip);
fprintf('  Motor: G(s) = %.1f/(%.1fs+1), sat 0-%d RPS\n', ...
        P.motor.K, P.motor.tau, P.motor.rps_max);
fprintf('  Hover: Omega0 = %.1f rad/s (%.0f RPS), T = %.2f N total\n', ...
        P.Omega0, P.Omega0/(2*pi), P.hover.T_total);