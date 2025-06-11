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
Nx = 200;
Ny = 200;
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
function f = make_frequency_axis(N, d)
if mod(N, 2) == 0
    f = (-N/2 : N/2 - 1) / (N * d);    % 중심 기준
else
    f = (-(N-1)/2 : (N-1)/2) / (N * d);
end
f = f * 2 * pi;   % rad/m 또는 rad/s
end

kx = make_frequency_axis(Nx, dx);
ky = make_frequency_axis(Ny, dy);
w  = make_frequency_axis(Nt, dt);

[Ky, Kx, W] = meshgrid(kx, ky, w);
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
for i = 1 : nFile
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

    disp([nFile, i]);
end
toc;