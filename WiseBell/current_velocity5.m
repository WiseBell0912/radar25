clear; close all; clc;

%% Search
file_path = 'D:/png2019/10/';
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
w = -pi/dt : 2*pi/Lt : pi/dt;
w(65) = [];

[Kx, Ky, W] = meshgrid(ky, kx, w);
K = sqrt(Kx.^2 + Ky.^2);

%% Windowing Function
window_xy = hann(Nx) * hann(Ny)';
window_t = hann(Nt);
window = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% Save
Date = NaT(length(file_list), 1);
SNR = zeros(length(file_list), 1);

%% Start parallel loop
tic
parfor i = 1 : length(file_list)

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

    %% Process

    % windowing
    image = png_wave;
    image_hanning = image .* window;

    % 3d-fft
    image_fft_SNR = fftn(image);
    image_fft_SNR = fftshift(image_fft_SNR);
    image_fft = fftn(image_hanning);
    image_fft = fftshift(image_fft);

    % image specturm
    image_spectrum_SNR = abs(image_fft_SNR).^2 / Nx^2 / Ny^2 / Nt^2;
    image_spectrum = abs(image_fft).^2 / Nx^2 / Ny^2 / Nt^2;

    % hpf
    K_th = 0.03;
    W_th = 0.35;
    hpf = (abs(W) > W_th & K > K_th);
    image_spectrum_hpf_SNR = image_spectrum_SNR .* hpf;
    image_spectrum_hpf = image_spectrum .* hpf;

    % weight
    E_th = 0.2 * max(image_spectrum_hpf, [], "all");
    wf = 2 * (image_spectrum_hpf > E_th & W > 0);
    weight = image_spectrum_hpf .* wf;

    % current surface
    a11 = sum(weight .* Kx .* Kx, "all");
    a12 = sum(weight .* Kx .* Ky, "all");
    a21 = sum(weight .* Ky .* Kx, "all");
    a22 = sum(weight .* Ky .* Ky, "all");
    A = [a11, a12; a21, a22];
    sigma = sqrt(g .* K .* tanh(K .* h));
    b11 = sum(weight .* Kx .* (W - sigma), "all");
    b21 = sum(weight .* Ky .* (W - sigma), "all");
    B = [b11; b21];
    U = A\B;
    Ux = U(1);
    Uy = U(2);

    % bpf
    dispersion_relationship = sigma + Kx.*Ux + Ky.*Uy;
    bpv1 = dispersion_relationship - 2 * 2*pi/Lt;
    bpv2 = dispersion_relationship + 2 * 2*pi/Lt;
    bpf = (bpv1 <= W & W <= bpv2);
    image_spectrum_bpf_SNR = image_spectrum_SNR .* bpf;
    image_spectrum_bpf = image_spectrum .* bpf;

    % snr
    signal = sum(image_spectrum_bpf_SNR, "all");
    noise = sum(image_spectrum_hpf_SNR, "all") - signal;
    snr = signal / noise;

    % result
    disp(['snr = ', num2str(snr)]);

    % save
    SNR(i) = snr;
    Date(i) = date;


    %% Checker
    disp([length(file_list), i]);

end
toc

save("snr10.mat", 'Date', 'SNR');