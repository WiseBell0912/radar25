clear; close all; clc;

%% 탐구 날짜 설정
load('date.mat');
date_find_val = datetime(2019, 10, 15, 11, 10, 00);
date_find_idx = find(Date == date_find_val);

%% Search
%file_path = '/mnt/usb/png2019/';
file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
file_list = dir([file_path, 'AIB_201910*.png']);

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

%% Input parameter
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
h_surf = 15;
h_wave = 30;

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
hpK_surf = (K > 0.0156810355); % surf
hpK_wave = (K > 0.0111906701); % wave

hpW = (abs(W) > 2*pi*0.03);

hpMask_surf = hpW .* hpK_surf;
hpMask_wave = hpW .* hpK_wave;

%%
nFile = length(file_list);

Date  = NaT(nFile,1);
surf_maxK = zeros(nFile, 1);
surf_maxW = zeros(nFile, 1);
wave_maxK = zeros(nFile, 1);
wave_maxW = zeros(nFile, 1);
surf_meanK = zeros(nFile, 1);
surf_meanW = zeros(nFile, 1);
wave_meanK = zeros(nFile, 1);
wave_meanW = zeros(nFile, 1);



for date_find_idx = 1 : 1

%% 파일 읽기 및 Zone 추출
png_path = fullfile(file_list(date_find_idx).folder, file_list(date_find_idx).name);
dateStr  = file_list(date_find_idx).name(5:end-4);
dateVal  = datetime(dateStr, 'InputFormat', 'yyyyMMdd_HHmm');

png_long = imread(png_path);
png_long = png_long(1 : 512*1080*128);
png_long = reshape(png_long, 512, 1080, 128);
png_long = flip(png_long, 2);
png_long = flip(png_long, 3);

[png_surf, png_wave, png_energy] = f_extract_zone_new( ...
    png_long, modifiy_theta, ...
    surf_theta1, surf_theta2, surf_center, ...
    wave_theta1, wave_theta2, wave_center, ...
    energy_theta1, energy_theta2, energy_center);

img_surf = single(png_surf);
img_wave = single(png_wave);
img_energy = single(png_energy);

%% 이미지 변환
img_surf = rot90(flip(flip(img_surf), 2));
img_wave = rot90(flip(flip(img_wave), 2));
img_energy = rot90(flip(flip(img_energy), 2));

%%
img_surf_spectrum = abs(fftn(img_surf .* window)).^2 / Nt^2 / Nx^2 / Ny^2;
img_surf_spectrum = fftshift(img_surf_spectrum);
%img_surf_spectrum_save(:, :, date_find_idx) = img_surf_spectrum(:, :, date_find_idx);
img_surf_spectrum = img_surf_spectrum .* hpMask_surf;

img_wave_spectrum = abs(fftn(img_wave .* window)).^2 / Nt^2 / Nx^2 / Ny^2;
img_wave_spectrum = fftshift(img_wave_spectrum);
%img_wave_spectrum_save(:, :, date_find_idx) = img_wave_spectrum(:, :, date_find_idx);
img_wave_spectrum = img_wave_spectrum .* hpMask_wave;

figure(1);
set(gcf, 'WindowState', 'maximized');
drawnow;
tiledlayout(2, 3);
sgtitle(datestr(dateVal, 'yyyy-mm-dd HH:MM'));
colormap gray;

nexttile;
surf(img_surf(:, :, 1), 'EdgeAlpha', 0);
title('surf');
view(0, 90);
axis equal;
axis off;
xlim([1, 201]);
ylim([1, 201]);

nexttile;
surf(Kx(:, :, 1), Ky(:, :, 1), sum(img_surf_spectrum, 3), 'EdgeAlpha', 0);
title('surf spectrum K');
view(0, 90);
axis equal;
colorbar;
xlim([min(Kx, [], 'all'), max(Kx, [], 'all')]);
ylim([min(Ky, [], 'all'), max(Ky, [], 'all')]);
zlim([0, 10]);
clim([0, 1]);       % 값이 0~1일 때만 색 변화

nexttile;
plot(w, squeeze(sum(sum(img_surf_spectrum, 1), 2)));
title('surf spectrum W');
ylim([0, 1]);

nexttile;
surf(img_wave(:, :, 1), 'EdgeAlpha', 0);
title('wave');
view(0, 90);
axis equal;
axis off;
xlim([1, 201]);
ylim([1, 201]);

nexttile;
surf(Kx(:, :, 1), Ky(:, :, 1), sum(img_wave_spectrum, 3), 'EdgeAlpha', 0);
title('wave spectrum K');
view(0, 90);
axis equal;
colorbar;
xlim([min(Kx, [], 'all'), max(Kx, [], 'all')]);
ylim([min(Ky, [], 'all'), max(Ky, [], 'all')]);
zlim([0, 10]);
clim([0, 1]);       % 값이 0~1일 때만 색 변화

nexttile;
plot(w, squeeze(sum(sum(img_wave_spectrum, 1), 2)));
title('wave spectrum W');
ylim([0, 1]);

png_name = ['./test/', datestr(dateVal, 'yyyy_mm_dd_HH_MM'), '.png'];
saveas(gcf, png_name);

disp(date_find_idx);

Date(date_find_idx) = dateVal;
surf_maxK(date_find_idx) = max(sum(img_surf_spectrum, 3), [], 'all');
surf_maxW(date_find_idx) = max(sum(sum(img_surf_spectrum, 1), 2), [], 'all');
wave_maxK(date_find_idx) = max(sum(img_wave_spectrum, 3), [], 'all');
wave_maxW(date_find_idx) = max(sum(sum(img_wave_spectrum, 1), 2), [], 'all');
surf_meanK(date_find_idx) = mean(sum(img_surf_spectrum, 3), 'all');
surf_meanW(date_find_idx) = mean(sum(sum(img_surf_spectrum, 1), 2), 'all');
wave_meanK(date_find_idx) = mean(sum(img_wave_spectrum, 3), 'all');
wave_meanW(date_find_idx) = mean(sum(sum(img_wave_spectrum, 1), 2), 'all');
end