syms phi theta psi;

% Вісь X - направлена вздовж носа - хвоста
% Вісь Y - направлена вгору
% Вісь Z - доповнює цю систему та спрямована в правий бік

% Перехід віз нормальної до зв'язаної

% Поворот навколо осі Zg нормальної системи координат тобто тангаж 
R_Zg_theta = [cos(theta) sin(theta)  0 ;
    -sin(theta) cos(theta)  0 ;
    0          0       1];

% Поворот навколо осі Yg нормальної системи координат тобто рискання 
R_Yg_psi = [cos(psi) 0 -sin(psi);
    0       1     0    ;
    sin(psi) 0 cos(psi)];

% Поворот навколо осі Xg нормальної системи координат тобто крен
R_Xg_phi = [1    0         0     ;
    0 cos(phi)  sin(phi) ;
    0 -sin(phi) cos(phi)];

% Повноцінна матриця повороту від нормальної системи координат до зв'язаної

R_normal2body = R_Xg_phi * R_Zg_theta * R_Yg_psi;

syms beta alpha;

% Перехід від швидкісної до зв'язаної

% Поворот навколо осі Ya швидкісної системи координат тобто кут ковзання

R_Ya_beta =   [cos(beta) 0 -sin(beta);
    0       1     0     ;
    sin(beta) 0 cos(beta)];

% Поворот навколо осі Za швидкісної системи координат тобто кут атаки

R_Za_alpha = [cos(alpha) sin(alpha)  0 ;
    -sin(alpha) cos(alpha)  0 ;
    0          0       1];

% Матриця переходу зі швидкісної до зв'язаної системи координат

R_speed2body = R_Za_alpha * R_Ya_beta;

syms Theta Psi;

% Перехід від нормальної до траєкторної

% Поворот навколо осі Yg нормальної системи координат на шляховий кут

R_Yg_Psi = [cos(Psi) 0 -sin(Psi);
    0      1     0    ;
    sin(Psi) 0 cos(Psi)];

% Поворот навколо осі Zg нормально системи координат на кут нахилу
% траєкторії

R_Zg_Theta = [cos(Theta) sin(Theta)  0 ;
    -sin(Theta) cos(Theta)  0 ;
    0          0       1];

% Матриця переходу від нормальної до траєкторної системи координат

R_normal2trajectory = R_Zg_Theta * R_Yg_Psi;

syms psi_a theta_a gamma_a;

% Матриця повороту від повітряної до швидкісної системи координат

R_air2speed = [1       0          0      ;
    0  cos(gamma_a)  sin(gamma_a) ;
    0 -sin(gamma_a)  cos(gamma_a)];

% Поворот навколо осі Yair повітряної системи координат на повітряний кут
% рискання

R_Yair_psi_a = [cos(psi_a) 0 -sin(psi_a);
    0      1     0      ;
    sin(psi_a) 0 cos(psi_a)];

% Поворот навколо осі Zair повітряної системи координат на повітряний кут
% тангажу

R_Zair_theta_a = [cos(theta_a) sin(theta_a)  0 ;
    -sin(theta_a) cos(theta_a)  0 ;
    0          0           1];

R_air2normal = R_Zair_theta_a * R_Yair_psi_a;

R_normal2air = R_air2normal';

R_normal2speed = R_air2speed * R_normal2air';
