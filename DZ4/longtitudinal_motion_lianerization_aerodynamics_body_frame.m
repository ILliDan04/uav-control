clear; clc;

%% ------------------------------------------------------------------------
%  1. SYMBOLIC VARIABLES
%  ------------------------------------------------------------------------

% --- States: body-frame linear velocities [m/s] --------------------------
syms U V W real

% --- States: body-frame angular rates [rad/s] ----------------------------
syms P Q R real

% --- States: Euler angles [rad] ------------------------------------------
syms Phi Theta real

% --- Controls: rotor angular speeds [rps] ------------------------------
syms w1 w2 w3 w4 real

% --- Mass / geometry / environment ---------------------------------------
syms m g rho S Lref b real      % Lref = reference length (L~), b = rotor arm

% --- Inertia tensor entries [kg*m^2] -------------------------------------
syms Ixx Iyy Izz Ixz real

% --- Aerodynamic coefficients --------------------------------------------
syms CA0 CAa2 CNa CMa CMq real  % CMq is the C_M^(omega~) damping derivative

% --- Propulsion ----------------------------------------------------------
syms CT d_prop real                    % thrust coefficient, d_prop = diameter of prop in m

% Assumptions: keep sqrt/atan from generating abs() and piecewise results
assumeAlso(m    > 0);
assumeAlso(rho  > 0);
assumeAlso(S    > 0);
assumeAlso(Lref > 0);
assumeAlso(Iyy  > 0);
assumeAlso(U    > 0);           % forward flight only


%% ------------------------------------------------------------------------
%  2. INTERMEDIATE (ALGEBRAIC) QUANTITIES
%  ------------------------------------------------------------------------

Vwind = sqrt(U^2 + V^2 + W^2);          % airspeed magnitude
alpha = atan(W/U);                      % angle of attack
qbar  = rho * Vwind^2 * S / 2;          % dynamic pressure


%% ------------------------------------------------------------------------
%  3. FORCE AND MOMENT BUILD-UP
%  ------------------------------------------------------------------------

% --- Aerodynamic drag and lift (wind-axis magnitudes) --------------------
Axial = (CA0 + CAa2*alpha^2) * qbar;
Normal = (CNa*alpha)          * qbar;

% --- Rotor thrust (all four rotors, axial) -------------------------------
Thrust = CT * rho * d_prop ^ 4 * (w1^2 + w2^2 + w3^2 + w4^2);

% --- Aerodynamic pitching moment (static + rate damping) -----------------
M_aero = (CMa*alpha + CMq*(Lref/Vwind)*Q) * qbar * Lref;

% --- Rotor differential-thrust pitching moment ---------------------------
M_rotor = CT * rho * d_prop ^ 4 * b * (w3^2 + w4^2 - w1^2 - w2^2);


%% ------------------------------------------------------------------------
%  4. EQUATIONS OF MOTION
%  ------------------------------------------------------------------------

% --- Axial velocity ------------------------------------------------------
Udot = -g*sin(Theta) - Q*W + R*V ...
       + (-Axial + Thrust ) / m;

% --- Normal velocity (no rotor contribution: rotors are purely axial) ----
Wdot =  g*cos(Phi)*cos(Theta) - P*V + Q*U ...
       - Normal / m;

% --- Pitch rate ----------------------------------------------------------
Qdot = ( P*R*(Izz - Ixx) + (R^2 - P^2)*Ixz + M_aero + M_rotor ) / Iyy;

% --- Pitch attitude ------------------------------------------------------
Thetadot = Q*cos(Phi) - R*sin(Phi);


%% ------------------------------------------------------------------------
%  5. ASSEMBLE VECTORS
%  ------------------------------------------------------------------------

f     = [ Udot ; Wdot ; Qdot ; Thetadot ];
x     = [ U    ; W    ; Q    ; Theta    ];
uctrl = [ w1   ; w2   ; w3   ; w4       ];
disp(f)

% Parameter list, for later subs() of numerical values
aoa_trim = 0.0548;
params = [ m g rho S Lref b Ixx Iyy Izz Ixz CA0 CAa2 CNa CMa CMq CT d_prop];
params_values = [2.7939 9.80665 1.225 0.0122718463 0.125 0.115 0.02 0.037 0.038 0.0 0.173 -3.56 8.200933208 -1.139274432 -5 0.0056 0.254];
initial_states = [U V W P Q R Phi Theta w1 w2 w3 w4];
initial_states_values_trim_aoa = [90*cos(aoa_trim) 0 90*sin(aoa_trim) 0 0 0 0 aoa_trim 252.3715 252.3715 369.1530 369.1530];
initial_states_values_0_aoa = [90 0 0 0 0 0 0 0 252.3715 252.3715 369.1530 369.1530];

wind_subs = subs(Vwind, initial_states, initial_states_values_trim_aoa);
alpha_subs = subs(alpha, initial_states, initial_states_values_trim_aoa);
qbar_subs = subs(qbar, [params, initial_states, Vwind], [params_values, initial_states_values_trim_aoa, wind_subs]);

disp(double(wind_subs));
disp(double(alpha_subs));
disp(double(qbar_subs));

%DERIVATIVES AS FUNCTIONS OF PARAMETERS

%Udot equation derivatives
UdotDu = subs(diff(Udot, U), initial_states, initial_states_values_0_aoa);
UdotDw = subs(diff(Udot, W), initial_states, initial_states_values_0_aoa);
UdotDq = subs(diff(Udot, Q), initial_states, initial_states_values_0_aoa);
UdotDtheta = subs(diff(Udot, Theta), initial_states, initial_states_values_0_aoa);
% disp("UdotDu: " + string(UdotDu));
% disp("UdotDw: " + string(UdotDw));
% disp("UdotDq: " + string(UdotDq));
% disp("UdotDtheta: " + string(UdotDtheta));

%Wdot equation derivatives
WdotDu = subs(diff(Wdot, U), initial_states, initial_states_values_0_aoa);
WdotDw = subs(diff(Wdot, W), initial_states, initial_states_values_0_aoa);
WdotDq = subs(diff(Wdot, Q), initial_states, initial_states_values_0_aoa);
WdotDtheta = subs(diff(Wdot, Theta), initial_states, initial_states_values_0_aoa);
% disp("WdotDu: " + string(WdotDu));
% disp("WdotDw: " + string(WdotDw));
% disp("WdotDq: " + string(WdotDq));
% disp("WdotDtheta: " + string(WdotDtheta));

%QDot equation derivatives
QdotDu = subs(diff(Qdot, U), initial_states, initial_states_values_0_aoa);
QdotDw = subs(diff(Qdot, W), initial_states, initial_states_values_0_aoa);
QdotDq = subs(diff(Qdot, Q), initial_states, initial_states_values_0_aoa);
QdotDtheta = subs(diff(Qdot, Theta), initial_states, initial_states_values_0_aoa);
% disp("QdotDu: " + string(QdotDu));
% disp("QdotDw: " + string(QdotDw));
% disp("QdotDq: " + string(QdotDq));
% disp("QdotDtheta: " + string(QdotDtheta));

%Thetadot equation derivatives
ThetadotDu = subs(diff(Thetadot, U), initial_states, initial_states_values_0_aoa);
ThetadotDw = subs(diff(Thetadot, W), initial_states, initial_states_values_0_aoa);
ThetadotDq = subs(diff(Thetadot, Q), initial_states, initial_states_values_0_aoa);
ThetadotDtheta = subs(diff(Thetadot, Theta), initial_states, initial_states_values_0_aoa);
% disp("ThetadotDu: " + string(ThetadotDu));
% disp("ThetadotDw: " + string(ThetadotDw));
% disp("ThetadotDq: " + string(ThetadotDq));
% disp("ThetadotDtheta: " + string(ThetadotDtheta));

disp("CHECK IF DYNAMICS ARE ZERO AT TRIM POINT");
fx0 = subs(f, params, params_values);
fx0 = subs(fx0, initial_states, initial_states_values_trim_aoa);
fx0 = double(fx0);
disp(fx0);

raw_A = jacobian(f, x);
raw_B = jacobian(f, uctrl);
Ax0 = subs(raw_A, initial_states, initial_states_values_trim_aoa);
Bx0 = subs(raw_B, initial_states, initial_states_values_trim_aoa);
Sim.A = double(subs(Ax0, params, params_values));
Sim.B = double(subs(Bx0, params, params_values));

Sim.C = eye(4);
Sim.D = zeros(4, 4);

Sim.trim_control = [252.3715 252.3715 369.1530 369.1530];

Sim.initial_conditions = [90*cos(aoa_trim) 90*sin(aoa_trim) 0.0 aoa_trim];

disp("A: ");
disp(Sim.A);
disp("B: ");
disp(Sim.B);
disp("Eig A: ");
disp(eig(Sim.A));