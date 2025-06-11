clear; close all; clc;
%% Search
% file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
file_path = '/media/zerog/EA168F94168F6085/png2020/';
file_list = dir([file_path, '*.png']);

%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(160);
surf_center = 870 + 600/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 600/2;

energy_theta1 = deg2rad(100);
energy_theta2 = deg2rad(185);
energy_center = 1750 + 500/2;

r = linspace(800, 2333, 512);
theta = linspace(0 - modifiy_theta, 2*pi - modifiy_theta, 1080);

x = r' * cos(theta);
y = r' * sin(theta);

%% input parameter
Nx = 201;
Ny = 201;
Nt = 128;

dx = 3;
dy = 3;
dt = 1.43;

Lx = 600;
Ly = 600;
Lt = Nt * dt;

g = 9.81;
surf_h = 30;
wave_h = 30;

%% Frequency
kx = -pi/dx : 2*pi/Lx : pi/dx;
ky = -pi/dy : 2*pi/Ly : pi/dy;
w  = -pi/dt : 2*pi/Lt : pi/dt;
w(65) = [];

[Kx, Ky, W] = meshgrid(ky, kx, w);
Kx = single(Kx);
Ky = single(Ky);
W  = single(W);
K  = sqrt(Kx.^2 + Ky.^2);

%% Windowing
window_xy = hann(Nx) * hann(Ny)';      % Nx x Ny
window_t  = hann(Nt);                  % Nt
window    = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% 하이패스 필터(HPF)용 마스크 미리 계산(반복문 밖)
hpK = (K > 0.03);
hpW = (abs(W) > 0.35);
hpMask = hpK & hpW;

%% 결과 저장용
nFile = length(file_list);

r_Date  = NaT(nFile,1);
r_LandE = zeros(nFile,1);
r_SurfE = zeros(nFile,1);
r_WaveE = zeros(nFile,1);

r_surf_K_max = zeros(nFile,1);
r_surf_K_mean = zeros(nFile,1);
r_surf_W_max = zeros(nFile,1);
r_surf_W_mean = zeros(nFile,1);

r_wave_K_max = zeros(nFile,1);
r_wave_K_mean = zeros(nFile,1);
r_wave_W_max = zeros(nFile,1);
r_wave_W_mean = zeros(nFile,1);

r_surf_SNR   = zeros(nFile,1);
r_surf_Ux    = zeros(nFile,1);
r_surf_Uy    = zeros(nFile,1);
r_surf_Tp   = zeros(nFile,1);
r_surf_Pdir = zeros(nFile,1);

r_wave_SNR   = zeros(nFile,1);
r_wave_Ux    = zeros(nFile,1);
r_wave_Uy    = zeros(nFile,1);
r_wave_Tp   = zeros(nFile,1);
r_wave_Pdir = zeros(nFile,1);

tic

parfor i = 1 : nFile
    %% 파일 읽기 및 Zone 추출
    png_path = fullfile(file_list(i).folder, file_list(i).name);
    dateStr  = file_list(i).name(5:end-4);
    dateVal  = datetime(dateStr, 'InputFormat', 'yyyyMMdd_HHmm');

    r_Date(i) = dateVal;

    png_long = imread(png_path);
    png_long = png_long(1 : 512*1080*128);
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    r_LandE(i) = sum(png_long(1:512, 355:385, :), 'all') + sum(png_long(1:512, 875:935, :), 'all');

    % figure(1);
    % clf;
    % axis equal;
    % hold on;
    % surf(x, y, png_long(:, :, 1), 'EdgeAlpha', 0);
    % surf(x(1:512, 355:385), y(1:512, 355:385), png_long(1:512, 355:385, 1), 'EdgeAlpha', 0.2);
    % surf(x(1:512, 875:935), y(1:512, 875:935), png_long(1:512, 875:935, 1), 'EdgeAlpha', 0.2);
    % hold off;

    disp(i);

end

save('landE.mat', 'r_Date', 'r_LandE');