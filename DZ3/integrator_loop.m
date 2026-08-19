clear L;

L.dt = 0.0001;
L.T = 10;

L.t = 0:dt:T;
L.N = numel(t);

L.z = zeros(1, N);

% [px, py, pz, vx, vy, vz, qw, qx, qy, qz, wx, wy, wz]
L.X0 = [0;0;0; 0;0;0; 1;0;0;0; 0;0;0;];
L.X = X0;


L.Fb = [0; 0; 0] .* (t >= 1);

L.Mb = [0;0;0;];

% Simulation loop
for k = 1:L.N
    L.X = integrator_step(L.X, Fb(:, k), L.Mb, L.dt, m, I);
    L.z(k) = L.X(3);
end

% Plot
figure;
plot(L.t, L.z)
xlabel('Time [s]')
ylabel('z')
grid on