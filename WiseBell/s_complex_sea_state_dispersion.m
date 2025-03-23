clear; close all; clc;

%% Information of Sea
Nx = 121;
Ny = 211;
Nt = 128;

dx = 3;
dy = 3;
dt = 1.43;

Lx = 360;
Ly = 630;
Lt = Nt * dt;

g = 9.81;  % 중력 가속도
h = 30;    % 수심 (깊은 물 가정)

%% Location
x = 0 : dx : Lx;
y = 0 : dy : Ly;
t = dt : dt : Lt;
[Y, X] = meshgrid(y, x);

%% Setup for multiple wave components
N_waves = 50;  % 복합 파동의 개수
w_min = 0;  % 최소 주파수
w_max = pi/dt;  % 최대 주파수

% 주파수, 각도, 초기 위상을 무작위로 생성
w = w_min + (w_max - w_min) * rand(N_waves, 1);  % 주파수 범위 내 무작위 값
theta = 2 * pi * rand(N_waves, 1);  % 각도 (0 ~ 2π 사이)
a = 1 * rand(N_waves, 1);               % 랜덤한 진폭
phi = -pi + 2 * pi * rand(N_waves, 1);  % 초기 위상 (-π ~ π 사이)

%% Sea State
n = zeros(Nx, Ny, Nt);
for n_wave = 1:N_waves
    % 파수 k 계산 (깊은 물 조건: tanh(kh) ≈ 1)
    k_wave = w(n_wave).^2 / g;

    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nt
                % 개별 파동을 계산하여 합산
                n(i, j, k) = n(i, j, k) + a(n_wave) * cos(k_wave * cos(theta(n_wave)) * x(i) + ...
                                                  k_wave * sin(theta(n_wave)) * y(j) - ...
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
