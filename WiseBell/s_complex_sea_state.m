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
h = 30;    % 수심

%% Location
x = 0 : dx : Lx;
y = 0 : dy : Ly;
t = dt : dt : Lt;
[Y, X] = meshgrid(y, x);

%% Setup: 여러 개의 파동
num_waves = 50;  % 합성할 파동의 개수

k_wave = 2 * pi * rand(1, num_waves) / 75;  % 랜덤한 파수
w = 2 * pi * rand(1, num_waves) / 100;       % 랜덤한 각 주파수
theta = 2 * pi * rand(1, num_waves);        % 랜덤한 파동 방향
a = 10 * rand(1, num_waves);               % 랜덤한 진폭
phi = 2 * pi * rand(1, num_waves) - pi;     % 랜덤한 초기 위상

%% Sea State: 여러 파동을 합성하여 복합 파동 생성
n = zeros(Nx, Ny, Nt);  % 해수면의 파동을 담을 배열
for wave = 1 : num_waves
    for i = 1 : Nx
        for j = 1 : Ny
            for k = 1 : Nt
                n(i, j, k) = n(i, j, k) + a(wave) * cos(k_wave(wave) * cos(theta(wave)) * x(i) + ...
                                                k_wave(wave) * sin(theta(wave)) * y(j) - w(wave) * t(k) + phi(wave));
            end
        end
    end
end

%% Display: 복합 파동 시각화
figure(1);
for i = 1 : Nt
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
