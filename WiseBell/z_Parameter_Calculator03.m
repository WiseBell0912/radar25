clear; close all; clc;

%% Search
%file_path = 'C:/Users/Hyeonjong Im/Documents/GitHub/Radar_2nd/PNG/';
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
Nx = 211;
Ny = 121;
Nt = 128;

dx = 3;
dy = 3;
dt = 1.43;

Lx = 630;
Ly = 360;
Lt = Nt * dt;

g = 9.81;

h = 30;

%% Windowing Function

% Method 1 (Hann)
window_xy = hann(Nx) * hann(Ny)';
window_t = hann(Nt);
window = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

% Method 2 (Tukey)
% r = 0.5;
% window_xy = tukeywin(Nx, r) * tukeywin(Ny, r)';
% window_t = tukeywin(Nt, r);
% window = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% Calculate Wave Parameter

% Set Storage
A_date = NaT(length(file_list), 1);
B_Tp_surf = zeros(length(file_list), 1);
C_Tp_wave = zeros(length(file_list), 1);
D_Pdir_surf = zeros(length(file_list), 1);
E_Pdir_wave = zeros(length(file_list), 1);
F_SNR_surf = zeros(length(file_list), 1);
G_SNR_wave = zeros(length(file_list), 1);
H_Land_energy = zeros(length(file_list), 1);

% Start
for i = 1 : length(file_list)
tic
    %% Radar Sub-Image sequence
    png_path = [file_list(i).folder '/' file_list(i).name];

    date = file_list(i).name(5:end-4);
    date = datetime(date, 'InputFormat', 'yyyyMMdd_HHmm');

    png_long = imread(png_path);
    png_long = png_long(1 : 512 * 1080 * 128);
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    [png_surf, png_wave] = f_extract_zone(png_long, modifiy_theta, surf_theta1, surf_theta2, surf_center, wave_theta1, wave_theta2, wave_center);

    %% Energy Level
    Land_energy = sum(png_long(:, 365:385, 1), 'all');

    %% Image Normalization
    png_surf = png_surf - mean(png_surf, 3);
    png_wave = png_wave - mean(png_wave, 3);

    % for j = 1:128
    %     s = mesh(png_wave(:, :, j), 'Edgecolor', 'flat');
    %     view(0, 90);
    %     rotate(s, [0, 0, 1], 90);
    %     axis equal;
    %     pause(0.1);
    % end

    %% Windowing
    png_surf = png_surf .* window;
    png_wave = png_wave .* window;

    clear var window_t window_xy

    %% FFT (with spectrum normalizaion)

    % For frequency
    image_spectra_surf_freq = abs(fftn(png_surf)).^2 ./ (Nx * Ny * Nt);
    image_spectra_wave_freq = abs(fftn(png_wave)).^2 ./ (Nx * Ny * Nt);

    image_spectra_surf_freq = fftshift(image_spectra_surf_freq);
    image_spectra_wave_freq = fftshift(image_spectra_wave_freq);

    image_spectra_surf_SNR = image_spectra_surf_freq;
    image_spectra_wave_SNR = image_spectra_wave_freq;

    % For physics
    image_spectra_surf_psic = abs(fftn(png_surf)).^2 ./ (Lx * Ly * Lt);
    image_spectra_wave_psic = abs(fftn(png_wave)).^2 ./ (Lx * Ly * Lt);

    image_spectra_surf_psic = fftshift(image_spectra_surf_psic);
    image_spectra_wave_psic = fftshift(image_spectra_wave_psic);

    %% Frequency
    kx = -pi/dx : 2*pi/Lx : pi/dx;
    ky = -pi/dy : 2*pi/Ly : pi/dy;
    w = -pi/dt : 2*pi/Lt : pi/dt;
    w(65) = [];

    [Kx, Ky, W] = meshgrid(ky, kx, w);

    Kx = single(Kx);
    Ky = single(Ky);
    W = single(W);

    K = sqrt(Kx.^2 + Ky.^2);

    %% HPF
    W_pass = 2 * pi * 0.03;
    filter_W = W_pass < abs(W);

    K_pass = 0.011191;
    filter_K = K_pass <= abs(K);

    filter_HP = filter_W .* filter_K;

    image_spectra_surf_freq = image_spectra_surf_freq .* filter_HP;
    image_spectra_wave_freq = image_spectra_wave_freq .* filter_HP;

    image_spectra_surf_psic = image_spectra_surf_psic .* filter_HP;
    image_spectra_wave_psic = image_spectra_wave_psic .* filter_HP;

    clear var W_pass K_pass filter_W filter_K filter_HP

    %% Current Velocity
    sigma = sqrt(g .* K .* tanh(K .* h));

    a1 = sum( image_spectra_surf_psic .* Kx.^2 , 'all');
    a2 = sum( image_spectra_surf_psic .* Ky.^2 , 'all');
    a3 = sum( image_spectra_surf_psic .* Kx .* Ky , 'all');
    a4 = sum( image_spectra_surf_psic .* (W - sigma) .* Kx.^2 , 'all');
    a5 = sum( image_spectra_surf_psic .* (W - sigma) .* Ky.^2 , 'all');

    b1 = sum( image_spectra_wave_psic .* Kx.^2 , 'all');
    b2 = sum( image_spectra_wave_psic .* Ky.^2 , 'all');
    b3 = sum( image_spectra_wave_psic .* Kx .* Ky , 'all');
    b4 = sum( image_spectra_wave_psic .* (W - sigma) .* Kx.^2 , 'all');
    b5 = sum( image_spectra_wave_psic .* (W - sigma) .* Ky.^2 , 'all');

    surf_ux = ( a2 * a4 - a3 * a5 ) / ( a1 * a2 - a3^2 );
    surf_uy = ( a1 * a5 - a3 * a4 ) / ( a1 * a2 - a3^2 );

    wave_ux = ( b2 * b4 - b3 * b5 ) / ( b1 * b2 - b3^2 );
    wave_uy = ( b1 * b5 - b3 * b4 ) / ( b1 * b2 - b3^2 );

   clear var a1 a2 a3 a4 a5 b1 b2 b3 b4 b5

   %% BPF
   filter_surf = abs(W - (sigma + Kx .* surf_ux + Ky .* surf_uy)) < 2 * (2*pi/Lt);
   filter_wave = abs(W - (sigma + Kx .* wave_ux + Ky .* wave_uy)) < 2 * (2*pi/Lt);

   image_spectra_surf_freq = image_spectra_surf_freq .* filter_surf;
   image_spectra_wave_freq = image_spectra_wave_freq .* filter_wave;

   image_spectra_surf_psic = image_spectra_surf_psic .* filter_surf;
   image_spectra_wave_psic = image_spectra_wave_psic .* filter_wave;

   clear var filter_surf filter_wave

   %% Modulation Transer Function
   image_spectra_surf_freq = image_spectra_surf_freq .* K.^(1.2);
   image_spectra_wave_freq = image_spectra_wave_freq .* K.^(1.2);

   image_spectra_surf_psic = image_spectra_surf_psic .* K.^(1.2);
   image_spectra_wave_psic = image_spectra_wave_psic .* K.^(1.2);

   %% Spectrum Analysis for Peak Period
   image_spectra_surf_1D = squeeze(sum(sum(image_spectra_surf_freq, 1), 2));
   image_spectra_wave_1D = squeeze(sum(sum(image_spectra_wave_freq, 1), 2));

   [~, idx_surf1] = max(image_spectra_surf_1D);
   [~, idx_wave1] = max(image_spectra_wave_1D);

   Tp_surf = 2*pi / abs(w(idx_surf1));
   Tp_wave = 2*pi / abs(w(idx_wave1));

   %% Spectrum Analysis for Wave Direction
   image_spectra_surf_2D = squeeze(sum(image_spectra_surf_psic, 3));
   image_spectra_wave_2D = squeeze(sum(image_spectra_wave_psic, 3));

   [~, idx_surf2] = max(image_spectra_surf_2D, [], "all");
   [~, idx_wave2] = max(image_spectra_wave_2D, [], "all");

   Pdir_surf = rad2deg(atan2(Ky(idx_surf2), Kx(idx_surf2)));
   Pdir_wave = rad2deg(atan2(Ky(idx_wave2), Kx(idx_wave2))); 

   %% Signal to Noise Ratio
   signal_surf = sum(image_spectra_surf_freq, 'all');
   signal_wave = sum(image_spectra_wave_freq, 'all');

   noise_surf = sum(image_spectra_surf_SNR, 'all') - signal_surf;
   noise_wave = sum(image_spectra_wave_SNR, 'all') - signal_wave;

   SNR_surf = signal_surf / noise_surf;
   SNR_wave = signal_wave / noise_wave;

   %% Save
   A_date(i) = date;
   B_Tp_surf(i) = Tp_surf;
   C_Tp_wave(i) = Tp_wave;
   D_Pdir_surf(i) = Pdir_surf;
   E_Pdir_wave(i) = Pdir_wave;
   F_SNR_surf(i) = SNR_surf;
   G_SNR_wave(i) = SNR_wave;
   H_Land_energy(i) = Land_energy;

   %% Checker
   disp([length(file_list), i]);
toc
end

save("parameter_ver3.mat", 'A_date', 'B_Tp_surf', 'C_Tp_wave', 'D_Pdir_surf', 'E_Pdir_wave', 'F_SNR_surf', 'G_SNR_wave');