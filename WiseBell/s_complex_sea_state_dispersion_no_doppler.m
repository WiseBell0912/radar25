clear; close all; clc;

%% Information of Sea
Nx = 201;
Ny = 201;
Nt = 128;

dx = 3;
dy = 3;
dt = 1.43;

Lx = 600;
Ly = 600;
Lt = Nt * dt;

g = 9.81;  % 중력 가속도
h = 50;    % 수심 (깊은 물 가정)

Ux = 3;
Uy = 0;

%% Location
x = 0 : dx : Lx;
y = 0 : dy : Ly;
t = dt : dt : Lt;
[Y, X] = meshgrid(y, x);

%% Setup for multiple wave components
N_waves = 50;  % 복합 파동의 개수
a_min = 1;   % 최소 진폭
a_max = 4;   % 최대 진폭
w_min = 0;     % 최소 주파수
w_max = pi/dt; % 최대 주파수

% 주파수, 각도, 초기 위상, 진폭을 무작위로 생성
rng(0);
w = 0.7 * ones(N_waves, 1);                             % 주파수 범위 내 무작위 값
theta = 2 * pi * rand(N_waves, 1);                  % 각도 (0 ~ 2π 사이)
k_wave = w.^2 / g;
kx = k_wave .* cos(theta);
ky = k_wave .* sin(theta);

phi = -pi + 2 * pi * rand(N_waves, 1);              % 초기 위상 (-π ~ π 사이)
a_n = a_min + (a_max - a_min) * rand(N_waves, 1);   % 진폭 범위 내 무작위 값

%% Sea State
n = zeros(Nx, Ny, Nt);
for n_wave = 1:N_waves
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nt
                % 개별 파동을 계산하여 합산 (진폭 a_n(n_wave)를 포함)
                n(i, j, k) = n(i, j, k) + a_n(n_wave) * cos(k_wave(n_wave) * cos(theta(n_wave)) * x(i) + ...
                                                       k_wave(n_wave) * sin(theta(n_wave)) * y(j) - ...
                                                       w(n_wave) * t(k) + phi(n_wave));
            end
        end
    end
end

%% Display
figure(1);
colormap(winter);  % 색상을 바다처럼 설정
for i = 1:Nt
    surf(X, Y, n(:, :, i), 'EdgeAlpha', 0);
    colormap('winter');
    zlim([-100, 100]);  % 파동의 크기를 보기 좋게 제한
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Height (m)');
    title(['Time = ', num2str(t(i), '%.2f'), ' s']);
    pbaspect([Lx Ly Lt]);
    colorbar;
    pause(0.1);  % 애니메이션 속도 조절
end
