function v_dot_body = velocity_dot_body(v_body, omega_body, total_force_body, m)
%   COMPUTE derivative of velocity in body frame
%   Use nonlinear equation to compute derivative of velocity in body
arguments (Input)
    v_body                  % Velocity in body frame
    omega_body              % Angular velocity in body frame
    total_force_body        % Forces in body frame
    m                       % mass
end

arguments (Output)
    v_dot_body      % Acceleration in body frame
end
    v_dot_body = total_force_body / m - cross(omega_body, v_body);
end