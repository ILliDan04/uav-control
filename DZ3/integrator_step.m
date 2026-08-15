function x1 = integrator_step(x0, Fb, Mb, h, m, J)
%   This function performs integration step with given initial x0, inputs
%   Fb, Mb, step size h and parameters m, J
%   x0 is structured as following:
%   px, py, pz     - positions in NED frame
%   vx, vy, vz     - velocities in BODY frame
%   qw, qx, qy, qz - quaternion from BODY to NED
%   wx, wy, wz     - angular velocities in BODY frame
%#codegen
arguments (Input)
    x0 (13, 1) double                   % Initial state [px, py, pz, vx, vy, vz, qw, qx, qy, qz, wx, wy, wz]
    Fb (3, 1)  double                   % Forces applied to body in body frame
    Mb (3, 1)  double                   % Torque applied to body in body frame
    h  (1, 1)  double {mustBeNonnegative}  % Step size
    m  (1, 1)  double {mustBePositive}  % Mass
    J  (3, 3)  double                   % Inertia matrix
end

arguments (Output)
    x1 (13, 1) double                   % Propagated state x1 = [px, py, pz, vx, vy, vz, qw, qx, qy, qz, wx, wy, wz]
                                        %   px, py, pz     - positions in NED frame
                                        %   vx, vy, vz     - velocities in BODY frame
                                        %   qw, qx, qy, qz - quaternion from BODY to NED
                                        %   wx, wy, wz     - angular velocities in BODY frame
end
    g  = [0; 0; 9.80665];               % NED, down positive
    p0 = x0(1:3);
    v0 = x0(4:6);
    q0 = x0(7:10);
    w0 = x0(11:13);

    % COMPUTE NEXT W USING RK4 FORM
    a1 = omega_dot_body(w0,                J, Mb);
    a2 = omega_dot_body(w0 + h * a1 * 0.5, J, Mb);
    a3 = omega_dot_body(w0 + h * a2 * 0.5, J, Mb);
    a4 = omega_dot_body(w0 + h * a3,       J, Mb);
    w1 = w0 + (h/6)*(a1 + 2*a2 + 2*a3 + a4);

    % MIDPOINT OMEGA USED TO INTEGRATE VELOCITY
    wm = w0 + (h/2)*a3;

    % PROPOGATE ATTITUDE USING EXPONENTIAL QUATERNION MAP
    dth = (h/2)*(w0 + w1) + (h^2/12)*cross(w0, w1);  % CONING COMPENSATION

    %MIDPOINT ATTITUDE CHANGE USED TO INTEGRATE VELOCITY LATER
    dthm = (h/4)*(w0 + w1) + (h^2/48)*cross(w0, w1);

    delta_q = quatexp(dth);
    delta_q_m = quatexp(dthm);

    q1 = quatmul(q0, delta_q);
    qm = quatmul(q0, delta_q_m);

    R0  = quat2dcm_local(q0);
    R1  = quat2dcm_local(q1);
    Rm  = quat2dcm_local(qm);

    % INTEGRATING VELOCITY WITH RK4
    % Total body forces = R * g * m + Fb, T - thrust, R - rotation from NED
    % to BODY
    total_forces_body_0 = R0' * m * g + Fb;
    total_forces_body_m = Rm' * m * g + Fb;
    total_forces_body_1 = R1' * m * g + Fb;

    b1 = velocity_dot_body(v0,                w0, total_forces_body_0, m);
    b2 = velocity_dot_body(v0 + h * b1 * 0.5, wm, total_forces_body_m, m);
    b3 = velocity_dot_body(v0 + h * b2 * 0.5, wm, total_forces_body_m, m);
    b4 = velocity_dot_body(v0 + h * b3,       w1, total_forces_body_1, m);

    v1 = v0 + (h/6)*(b1 + 2*b2 + 2*b3 + b4);

    % INTEGRATE POSITION WITH RK4
    p1 = p0 + (h/6)*(R0*v0 + 2*Rm*(v0 + h/2*b1) + 2*Rm*(v0 + h/2*b2) + R1*(v0 + h*b3));

    x1 = [p1; v1; q1; w1];
end