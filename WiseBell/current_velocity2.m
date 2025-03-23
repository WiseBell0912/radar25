clear; close all; clc;

%% Search
file_path = 'D:/png2019/09/';
file_list = dir([file_path, '*.png']);

%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(165);
surf_center = 830 + 630/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 630/2;

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

%% Start parallel loop
tic
for i = 3239:3239%1 : length(file_list)

    %% Radar Sub-Image sequence
    png_path = [file_list(i).folder '/' file_list(i).name];
    date = file_list(i).name(5:end-4);
    date = datetime(date, 'InputFormat', 'yyyyMMdd_HHmm');

    png_long = imread(png_path);
    png_long = png_long(1:512*1080*128);
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    [png_surf, png_wave] = f_extract_zone(png_long, modifiy_theta, surf_theta1, surf_theta2, surf_center, wave_theta1, wave_theta2, wave_center);

    %%
    % 주어진 파라미터
    Nx = 121;  % x 방향 픽셀 수
    Ny = 211;  % y 방향 픽셀 수
    Nt = 128;  % 시간 축 픽셀 수

    dx = 3;    % x 방향 픽셀 당 물리적 거리 (m)
    dy = 3;    % y 방향 픽셀 당 물리적 거리 (m)
    dt = 1.43; % 시간축 픽셀 당 시간 간격 (s)

    Lx = 360;  % 전체 x 방향 물리적 크기 (m)
    Ly = 630;  % 전체 y 방향 물리적 크기 (m)
    Lt = Nt * dt; % 전체 시간 길이 (s)

    g = 9.81;  % 중력 가속도 (m/s^2)
    h = 30;    % 수심 (m)

    % 레이더 이미지 데이터 (예시로 랜덤 데이터 사용)
    % 실제 데이터가 있으면 그 데이터를 여기에 로드해야 합니다.
    eta = png_wave; % x, y, t 도메인에서 레이더 이미지

    % 푸리에 변환 수행 (3차원 푸리에 변환)
    F_eta = fftn(eta); % 3차원 푸리에 변환

    % 주파수 및 파수 벡터 계산
    kx = 2 * pi * (-floor(Nx/2):ceil(Nx/2)-1) / (Nx * dx);  % x 방향 파수
    ky = 2 * pi * (-floor(Ny/2):ceil(Ny/2)-1) / (Ny * dy);  % y 방향 파수
    omega = 2 * pi * (-floor(Nt/2):ceil(Nt/2)-1) / (Nt * dt); % 시간 방향 주파수

    % meshgrid로 kx와 ky를 2차원 격자로 변환
    [kx_grid, ky_grid] = meshgrid(kx, ky);

    % 푸리에 변환된 데이터 중앙으로 정렬
    F_eta_shifted = fftshift(F_eta);

    % 크기 정보
    k_magnitude = sqrt(kx_grid.^2 + ky_grid.^2);  % 파수의 크기 계산

    % 결과 시각화
    figure;
    subplot(2,2,1);
    imagesc(kx, ky, abs(squeeze(sum(F_eta_shifted, 3)))); % 시간축에 대해 합산하여 2D 플롯
    title('공간 주파수 도메인 (k_x, k_y)');
    xlabel('k_x (1/m)');
    ylabel('k_y (1/m)');
    colorbar;

    subplot(2,2,2);
    plot(omega, abs(squeeze(sum(sum(F_eta_shifted, 1), 2)))); % 공간축에 대해 합산하여 1D 플롯
    title('시간 주파수 도메인 (\omega)');
    xlabel('\omega (1/s)');
    ylabel('Amplitude');
    grid on;


    %% Checker
    disp([length(file_list), i]);

end
toc