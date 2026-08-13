function w_dot_body = omega_dot_body(w_body, J_body, M_body)
%   COMPUTES angular acceleration in body frame 
%   computes angular acceleration using nonlinear equation of rotational
%   motion
arguments (Input)
    w_body               % Angular velocity in body frame
    J_body               % Inertia matrix
    M_body               % Applied torque in body frame
end

arguments (Output)
    w_dot_body           % Output angular acceleration in body frame
end
    rhs = (M_body - cross(w_body, J_body*w_body));
    w_dot_body = J_body \ rhs;
end