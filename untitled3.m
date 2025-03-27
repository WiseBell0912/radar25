%% Sea Surface Video Simulation using JONSWAP Spectrum with Uniform Current
% 이 코드는 JONSWAP 스펙트럼 기반 해수면을 동적으로 생성하고,
% 균일 해류 효과를 포함하여 해수면의 시간 변화(동영상)를 저장하며,
% 동시에 각 시간 단계의 해수면 데이터(eta)를 MAT 파일로 저장합니다.

clear; close all; clc;

%% 1. 도메인 및 격자 설정
Nx = 201;         % x 방향 격자 포인트 수
Ny = 201;         % y 방향 격자 포인트 수
Lx = 600;         % x 방향 도메인 길이 (m)
Ly = 600;         % y 방향 도메인 길이 (m)

dx = 3;
dy = 3;
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

%% 2. 파수 그리드 생성
dkx = 2*pi/Lx;
dky = 2*pi/Ly;
kx = (-Nx/2:Nx/2-1)*dkx;
ky = (-Ny/2:Ny/2-1)*dky;
[KX, KY] = meshgrid(kx, ky);
K = sqrt(KX.^2 + KY.^2);

%% 3. JONSWAP 스펙트럼 파라미터 및 계산 (깊은 바다 가정)
g = 9.81;         % 중력가속도 (m/s^2)
Tp = 8;           % 피크 주기 (s)
omega_p = 2*pi/Tp; % 피크 각주파수 (rad/s)
Hs = 2;           % 유의 파고 (m)
gamma = 3.3;      % 피크 강화 인자
alpha = 0.076 * (Hs^2 * omega_p^4 / g^2);  % 스펙트럼 스케일링 파라미터

% 깊은 바다에서 dispersion relation: omega = sqrt(g*k)
omega = sqrt(g*K);
omega(K==0) = 0;  % k = 0 성분은 0으로 처리

% sigma: omega <= omega_p이면 0.07, 초과이면 0.09
sigma_val = zeros(size(K));
idx = (K > 0);
sigma_val(idx) = (omega(idx) <= omega_p)*0.07 + (omega(idx) > omega_p)*0.09;

% JONSWAP 스펙트럼 S(omega)
S_omega = zeros(size(K));
S_omega(idx) = (alpha*g^2) ./ (omega(idx).^5) .* exp(-5/4*(omega_p./omega(idx)).^4) .* ...
    (gamma.^( exp(-((omega(idx)-omega_p).^2)./(2*(sigma_val(idx)*omega_p).^2) )));
  
% omega와 k의 미분: domega/dk = 0.5*sqrt(g./K)
S_k = zeros(size(K));
S_k(idx) = S_omega(idx) .* (0.5 * sqrt(g./K(idx)));
S_k(K==0) = 0;

%% 4. 균일 해류 파라미터 (Uniform Current)
Ux = 10;   % x 방향 해류 속도 (m/s)
Uy = 0;    % y 방향 해류 속도 (m/s)

%% 5. 초기 복소 스펙트럼 생성 (랜덤 위상 부여)
Amplitude = sqrt(S_k);         % 각 파수 성분의 진폭
phi = 2*pi * rand(Ny, Nx);       % 0 ~ 2pi 사이의 랜덤 위상
A_k = Amplitude .* exp(1i*phi);  % 복소 스펙트럼 생성

%% 6. 시간 파라미터 및 동적 시뮬레이션
T_total = 183.04;      % 전체 시뮬레이션 시간 (s)
dt = 1.43;             % 시간 간격 (s)
Nt = 128;
time_vec = (0:Nt-1)*dt;

% 해수면 데이터를 저장할 3D 배열 (Ny x Nx x Nt)
eta_all = zeros(Ny, Nx, Nt);

%% 7. Video Writer 설정
v = VideoWriter('sea_surface_simulation.mp4', 'MPEG-4');
v.FrameRate = 10;  % 초당 프레임 수
open(v);

figure('Name','Sea Surface Video','NumberTitle','off');
for t = 1:Nt
    current_time = time_vec(t);
    % 시간 진화 인자: 해양 파동의 분산관계와 균일 해류 효과 반영
    time_factor = exp(-1i * (omega + KX*Ux + KY*Uy) * current_time);
    A_k_t = A_k .* time_factor;

    % 역 FFT를 통해 공간 영역 해수면 계산
    eta = ifft2(ifftshift(A_k_t), 'symmetric');

    % 해수면의 RMS 정규화 (일반적으로 Hs ≈ 4*RMS)
    rms_current = std(eta(:));
    desired_rms = Hs / 4;
    eta = eta * (desired_rms / rms_current);
    
    % 생성된 해수면 데이터를 저장 (각 시간 단계)
    eta_all(:,:,t) = eta;

    % 해수면 플롯 (영상용)
    surf(X, Y, eta, 'EdgeColor','none');
    shading interp;
    colormap winter;
    xlabel('x (m)'); ylabel('y (m)'); zlabel('Elevation (m)');
    title(sprintf('Time = %.1f s, U = [%.2f, %.2f] m/s', current_time, Ux, Uy));
    axis([0 Lx 0 Ly -Hs/2 Hs/2]);
    view(0, 90);
    colorbar;
    drawnow;
    pause(0.1);

    % 프레임 캡처 후 영상에 기록
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Video saved as sea_surface_simulation.mp4\n');

%% 8. 해수면 데이터 저장
% 예: eta_all, time_vec, x, y 등의 변수 저장
save('sea_surface_data.mat', 'eta_all', 'time_vec', 'x', 'y');
fprintf('Sea surface data saved as sea_surface_data.mat\n');
