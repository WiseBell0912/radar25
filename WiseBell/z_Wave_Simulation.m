clear; clc; close all;

%% Parameters
Nx = 50; % x 방향의 파동 수 (더 촘촘하게 설정)
Ny = 50; % y 방향의 파동 수 (더 촘촘하게 설정)
Lx = 10; % x 방향의 도메인 크기 [m]
Ly = 10; % y 방향의 도메인 크기 [m]
g = 9.81; % 중력 가속도 [m/s^2]
h = 30; % 수심 [m]
U = [1, 0]; % 해류 속도 [Ux, Uy] [m/s]
t = 0:1:100; % 시간 [s], 더 촘촘하게 설정

% x, y 도메인 설정
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% 파동 파라미터 설정
A = rand(Nx, Ny) * 0.5; % 각 파동의 진폭 (임의의 값)
phi = 2 * pi * rand(Nx, Ny); % 각 파동의 무작위 위상 [rad]

% 파수 벡터 계산
kx = 2 * pi * (0:Nx-1) / Lx; % x 방향 파수 벡터
ky = 2 * pi * (0:Ny-1) / Ly; % y 방향 파수 벡터

% 그리드 상에서 k 벡터 생성
[Kx, Ky] = meshgrid(kx, ky);
K = sqrt(Kx.^2 + Ky.^2); % 파수 벡터의 크기

% 주파수 계산 (분산 관계 사용)
omega = sqrt(g * K .* tanh(K * h)) + Kx * U(1) + Ky * U(2);

%% Wave Elevation Calculation
% 시뮬레이션을 위한 파동 계산
eta_total = zeros(size(X, 1), size(Y, 2), length(t)); % 파동의 총 합

for ti = 1:length(t)
    ti
    eta = zeros(size(X)); % 각 시간에서의 파동
    
    % 각 파수에 대해 파동을 더함
    for m = 1:Nx
        for n = 1:Ny
            % 파동 계산
            k_dot_r = Kx(m, n) * X + Ky(m, n) * Y; % k · r
            omega_t = omega(m, n) * t(ti); % ω * t
            eta = eta + A(m, n) * cos(k_dot_r - omega_t - phi(m, n)); % η 계산
        end
    end
    
    % 시간에 따른 해수면 변위 저장
    eta_total(:, :, ti) = eta;
end

%% Plotting the Result
figure;
for ti = 1:length(t)
    surf(X, Y, eta_total(:, :, ti));
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Wave Elevation [m]');
    title(['Sea Surface Elevation at t = ', num2str(t(ti)), ' s']);
    axis([0 Lx 0 Ly -2 2]);
    shading interp;
    view(2); % 위에서 본 시각으로 설정
    colorbar;
    pause(0.05); % 애니메이션 효과를 위한 일시 정지
end
