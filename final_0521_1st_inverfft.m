clear; close all; clc;

%% Search
%file_path = 'E:/png2019/10/';
file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
%file_path = '/media/zerog/EA168F94168F6085/png2019/';

%file_list = dir([file_path, '*.png']);
file_list = dir(fullfile(file_path, '*_201910*.png'));

%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(160);
surf_center = 870 + 600/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 600/2;

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
h = 30;

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

%% Save
nFile = length(file_list);

r_Date  = NaT(nFile,1);

r_LandE   = zeros(nFile,1);

r_surf_Singal = zeros(nFile,1);
r_wave_Signal = zeros(nFile,1);
r_surf_Noise  = zeros(nFile,1);
r_wave_Noise  = zeros(nFile,1);

%% Loop
tic
for i = 696
    % ---------- Read png ---------- %
    png_path = fullfile(file_list(i).folder, file_list(i).name);
    dateStr  = file_list(i).name(5:end-4);
    Date_Val  = datetime(dateStr, 'InputFormat', 'yyyyMMdd_HHmm');

    png_long = imread(png_path);
    png_long = png_long(1 : 512*1080*128);
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    [png_surf, png_wave] = f_extract_zone( ...
        png_long, modifiy_theta, ...
        surf_theta1, surf_theta2, surf_center, ...
        wave_theta1, wave_theta2, wave_center);

    % ---------- Edit image ---------- %
    % surf
    img_surf = single(png_surf);
    img_surf = rot90(flip(img_surf));
    % wave
    img_wave = single(png_wave);
    img_wave = rot90(flip(img_wave));

    % ---------- Land energy ---------- %
    LandE_Val = sum(png_long(1:512, 355:385, :), 'all') + sum(png_long(1:512, 875:935, :), 'all');

    % ---------- 3D FFT ---------- %
    % surf
    spectrum_surf_normal_non = fftshift(fftn(img_surf));
    % wave
    spectrum_wave_normal_non = fftshift(fftn(img_wave));

    % ---------- Spectrum ---------- %
    % surf
    spectrum_surf_power_non = abs(spectrum_surf_normal_non).^2;
    % wave
    spectrum_wave_power_non = abs(spectrum_wave_normal_non).^2;

    % ---------- HPF ---------- %
    hpK = (0.0231642599 <= K) & (K <= 0.1609926865);
    hpW = (2*pi/17 <= abs(W)) & (abs(W) <= 2*pi/5);
    hpMask = hpK & hpW;
    % surf
    spectrum_surf_power_non_hp = spectrum_surf_power_non .* hpMask;
    % wave
    spectrum_wave_power_non_hp = spectrum_wave_power_non .* hpMask;

    % ---------- SNR ---------- %
    % surf
    surf_signal = sum( spectrum_surf_power_non_hp(:) );
    surf_noise  = sum( spectrum_surf_power_non(:) ) - surf_signal;
    % wave
    wave_signal = sum( spectrum_wave_power_non_hp(:) );
    wave_noise  = sum( spectrum_wave_power_non(:) ) - wave_signal;

    % ---------- Save ---------- %
    r_Date(i)  = Date_Val;

    r_LandE(i) = LandE_Val;

    r_surf_Singal(i) = surf_signal;
    r_wave_Signal(i) = wave_signal;
    r_surf_Noise(i) = surf_noise;
    r_wave_Noise(i) = wave_noise;

    % ---------- Inverse FFT ---------- %
    fft = fftn(img_surf - mean(img_surf, 'all'));
    spectrum1 = abs(fft).^2;
    [~, idx_max] = max(spectrum1, [], "all");
    [idx_max_ky, idx_max_kx, idx_max_w] = ind2sub(size(spectrum1), idx_max);

    figure(1);
    tiledlayout(2, 1);
    nexttile;
    plot(squeeze(spectrum1(idx_max_ky, idx_max_kx, :)));
    xlabel('w');
    nexttile;
    plot(squeeze(spectrum1(:, idx_max_kx, idx_max_w)));
    xlabel('kx');


    fft = fftshift(fft);
    spectrum2 = abs(fft).^2;
    [~, idx_max] = max(spectrum2, [], "all");
    [idx_max_ky, idx_max_kx, idx_max_w] = ind2sub(size(spectrum2), idx_max);

    figure(2);
    tiledlayout(2, 1);
    nexttile;
    plot(squeeze(spectrum2(idx_max_ky, idx_max_kx, :)));
    xlabel('w');
    nexttile;
    plot(squeeze(spectrum2(:, idx_max_kx, idx_max_w)));
    xlabel('kx');

    % fft_power = abs(fft).^2;
    % fft = fft .* hpMask;
    % fft_power_hp = abs(fft).^2;
    % fft = ifftshift(fft);
    % img_surf_inverse = real(ifftn(fft));
    %
    % [~, idx_max] = max(fft_power_hp, [], "all");
    % [idx_max_ky, idx_max_kx, idx_max_w] = ind2sub(size(fft_power), idx_max);
    %
    % spectrum = squeeze(fft_power_hp(:, idx_max_ky, :));
    %
    % for j = 1 : 128
    %
    %     figure(1);
    %
    %     tiledlayout(2, 3);
    %     nexttile;
    %     imagesc(img_surf(:, :, j));
    %     nexttile;
    %     plot(w, log(squeeze(fft_power(idx_max_ky, idx_max_kx, :))));
    %     xlabel('w');
    %     xlim([min(w), max(w)]);
    %     nexttile;
    %     plot(kx, log(squeeze(fft_power(:, idx_max_kx, idx_max_w))));
    %     xlabel('kx');
    %     xlim([min(kx), max(kx)]);
    %
    %     nexttile;
    %     imagesc(img_surf_inverse(:, :, j));
    %     nexttile;
    %     plot(w, log(squeeze(fft_power_hp(idx_max_ky, idx_max_kx, :))));
    %     xlabel('w');
    %     xlim([min(w), max(w)]);
    %     nexttile;
    %     plot(kx, log(squeeze(fft_power_hp(:, idx_max_kx, idx_max_w))));
    %     xlabel('kx');
    %     xlim([min(kx), max(kx)]);
    %
    %     pause(0.1);
    %
    % end

    disp([nFile, i]);
end
toc;

%% Save
% save( ...
%     "snr_y2019_12_ver0520.mat", ...
%     'r_Date', ...
%     'r_LandE', ...
%     'r_surf_Singal', ...
%     'r_wave_Signal', ...
%     'r_surf_Noise', ...
%     'r_wave_Noise'...
%     );