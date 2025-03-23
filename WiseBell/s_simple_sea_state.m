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

g = 9.81;

h = 30;

%% Location
x = 0 : dx : Lx;
y = 0 : dy : Ly;
t = dt : dt : Lt;
[Y, X] = meshgrid(y, x);

%% Setup
k_wave = 2 * pi / 75;  % 파수
w = 2 * pi / 100;       % 각 주파수
theta = pi / 6;        % 파동 방향
a = 10;               % 진폭
phi = 0;               % 초기 위상

%% Sea State
n = zeros(Nx, Ny, Nt);
for i = 1 : Nx
    for j = 1 : Ny
        for k = 1 : Nt
            n(i, j, k) = a * cos(k_wave * cos(theta) * x(i) + k_wave * sin(theta) * y(j) - w * t(k) + phi);
        end
    end
end

%% Display
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
