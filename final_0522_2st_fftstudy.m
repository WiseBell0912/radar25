clear; close all; clc;

%% 1. 크기 및 분포 설정
Nx = 64; Ny = 64; Nt = 64;     % 크기
sigma = 1;                     % Rayleigh scale parameter

%% 2. 진폭 (Rayleigh), 위상 (Uniform)
amp_half = raylrnd(sigma, [Nx, Ny, Nt/2 + 1]);        % 진폭: 단측 스펙트럼
phi_half = 2 * pi * rand(Nx, Ny, Nt/2 + 1);           % 위상: 0~2pi

F_half = amp_half .* exp(1i * phi_half);              % 복소수 스펙트럼

%% 3. Hermitian symmetry 구성 (ifftn 위해 필요)
F_full = zeros(Nx, Ny, Nt, 'like', F_half);           % 전체 스펙트럼

% 절반 복사
F_full(:,:,1:Nt/2+1) = F_half;

% 공액 대칭 생성: z축 기준 반대쪽 복사 (2~Nt/2 → Nt:-1:Nt/2+2)
for k = 2:Nt/2
    F_full(:,:,Nt - k + 2) = conj(F_half(:,:,k));
end

F_full(:,:,1) = 0;               % DC 성분 제거
F_full(:,:,Nt/2+1) = 0;          % Nyquist 주파수도 실수 또는 제거

%% 4. 3D 역푸리에 변환으로 시공간 해수면 생성
eta = real(ifftn(F_full));       % 실수 시공간 파동장

%% 5. Parseval 정리 검증
energy_space = sum(abs(eta(:)).^2) / (Nx*Ny*Nt);
energy_freq  = sum(abs(F_full(:)).^2) / (Nx*Ny*Nt)^2;

fprintf('3D Parseval check:\n');
fprintf(' - Space domain energy : %.6f\n', energy_space);
fprintf(' - Freq domain energy  : %.6f\n', energy_freq);
fprintf(' - Difference           : %.2e\n', abs(energy_space - energy_freq));

%% 6. 시각화 (t = 1 시점)
figure;
for i = 1 : Nt
surf(eta(:,:,i), 'EdgeAlpha', 0);
title('\eta(x,y,t=1) from 3D Rayleigh Spectrum');
xlabel('x'); ylabel('y');
zlim([0, 0.2])
pause(0.5)
end