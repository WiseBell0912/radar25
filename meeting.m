clear;  clc;

%% ---------- 기본 파라미터 ---------- %%
g = 9.81;
h = 30;
dx = 3; dy = 3; dt = 1.43;

% 원 신호 크기
Nx = 201; Ny = 201; Nt = 128;
Lx = Nx * dx;
Ly = Ny * dy;
Lt = Nt * dt;

% 패딩 크기 (2의 거듭제곱 추천)
Nxp = 256; Nyp = 256; Ntp = 256;

%% ---------- 각파수 축 (FFT 기준, rad/m, rad/s) ---------- %%
dkx = 2*pi / (Nxp * dx);
dky = 2*pi / (Nyp * dy);
dw  = 2*pi / (Ntp * dt);

kx = (-floor(Nxp/2) : ceil(Nxp/2)-1) * dkx;
ky = (-floor(Nyp/2) : ceil(Nyp/2)-1) * dky;
w  = (-floor(Ntp/2) : ceil(Ntp/2)-1) * dw;

[Kx, Ky, W] = meshgrid(kx, ky, w);  % 주의: ky, kx, w 순서
K = sqrt(Kx.^2 + Ky.^2);

%% ---------- 단일파 생성 ---------- %%
% 원래 신호 그리드
x = (0:Nx-1) * dx;
y = (0:Ny-1) * dy;
t = (0:Nt-1) * dt;
[X, Y, T] = meshgrid(x, y, t);  % y, x, t 순서

% 파라미터
kx1 = 2 * pi / 100;      % rad/m
ky1 = 0;
K1 = sqrt(kx1^2 + ky1^2);
ux1 = 14;               % m/s
uy1 = 0;
w1 = sqrt(g * K1 * tanh(K1 * h)) + kx1 * ux1;
phi = pi;               % rad

% 신호
n = 3/2 * cos(kx1 * X + w1 * T) + 2 * cos(kx1 * X + w1 * T + phi);

figure;
imagesc(n(:, :, 1));

%% ---------- Windowing ---------- %%
window_xy = hann(Ny) * hann(Nx)';
window_t  = hann(Nt);
window = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% ---------- Padding ---------- %%
n_pad = zeros(Nyp, Nxp, Ntp);
n_pad(1:Ny, 1:Nx, 1:Nt) = n;

window_pad = zeros(Nyp, Nxp, Ntp);
window_pad(1:Ny, 1:Nx, 1:Nt) = window;

n_win = n_pad;% .* window_pad;

%% ---------- 고역통과 필터 ---------- %%
hpK = (2*pi/15 <= K) & (K <= 2*pi/5);
hpW = (2*pi/30 <= abs(W));
hpMask = hpK & hpW;

%% ---------- FFT & 스펙트럼 ---------- %%
spectrum = (fftn(n_win - mean(n_win, 3)));
spectrum = abs(spectrum).^2 / 256^2 / 256^2 / 256^2;% .* hpMask;

%% ---------- Threshold 적용 ---------- %%
threshold = 0 * max(spectrum(:));
filter = (spectrum >= threshold);

%% ---------- 분산관계 계산 ---------- %%
sigma = sqrt(g * K .* tanh(K * h));

%% ---------- Least Squares 유속 추정 ---------- %%
a11 = sum(Kx(filter).^2);
a12 = sum(Kx(filter) .* Ky(filter));
a22 = sum(Ky(filter).^2);
b1 = sum((W(filter) - sigma(filter)) .* Kx(filter));
b2 = sum((W(filter) - sigma(filter)) .* Ky(filter));

det = a11*a22 - a12^2;
ux = ( a22 * b1 - a12 * b2 ) / det;
uy = ( a11 * b2 - a12 * b1 ) / det;

fprintf('\n=== 유속 추정 결과 ===\n');
fprintf('Ux: %.3f m/s\n', ux);
fprintf('Uy: %.3f m/s\n', uy);
fprintf('속도 크기: %.3f m/s\n', sqrt(ux^2 + uy^2));

%% ---------- 시각화 ---------- %%
figure;
surf(n(:, :, 1));

% ω 방향 스펙트럼
figure;
plot(w, squeeze(sum(sum(spectrum, 1), 2)));
xlabel('\omega [rad/s]'); ylabel('Spectral Energy');
title('\omega-directional Spectrum');

% kx-ky 2D 스펙트럼
[~, idx_max] = max(spectrum(:));
[ky_idx, kx_idx, w_idx] = ind2sub(size(spectrum), idx_max);

figure;
plot(spectrum(:, kx_idx, w_idx));
view(0, 90); axis tight;
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');
title('2D Spectrum at Peak k');
colormap(jet); colorbar;