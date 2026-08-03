clear;
clc;

params.M = 0.5;
params.m = 0.2;
params.b = 0.1;
params.I = 0.006;
params.g = 9.8;
params.l = 0.3;

% Step force input
params.U = @(t) 1*((t >= 0.2) & (t < 0.25));

% -- SIMULATION -- 

% Initial state
X0 = [0; 0; 0; 0];

%tspan = 0:0.01:10;
tspan = [0,1];
[t,X] = ode45(@(t,X) invertedPendulumODE(t,X,params), tspan, X0);

% Plot cart position
figure;
plot(t,X(:,1),'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Cart position x (m)');
grid on;


% Plot pendulum angle
figure;
plot(t,X(:,3),'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Pendulum angle \phi (rad)');
grid on;

% ---------
% FUNCTIONS
% ---------

function dX = invertedPendulumODE(t,X,params)
% Params
M = params.M;
m = params.m;
b = params.b;
I = params.I;
g = params.g;
l = params.l;

u = params.U(t);

% State
% x      = X(1); Unused
xdot   = X(2);
phi    = X(3);
phidot = X(4);

A = [
    -m*l, I+m*l^2;
    M+m,  -m*l;
    ];

B = [
    m*g*l*phi;
    u-b*xdot;
    ];

q = A\B;

xddot = q(1);
phiddot = q(2);

% state derivative
dX = [
    xdot;
    xddot;
    phidot;
    phiddot;
];

end