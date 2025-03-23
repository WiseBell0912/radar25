clear; close all; clc;

%% Search
file_path = [pwd, '/PNG/'];
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

Lx = 600;
Ly = 600;
Lt = Nt * dt;

g = 9.81;

h = 30;

%% Frequency
kx = -pi/dx : 2*pi/Lx : pi/dx;
ky = -pi/dy : 2*pi/Ly : pi/dy;
w = -pi/dt : 2*pi/Lt : pi/dt;
w = fftshift(w);
w(65) = [];

[Kx, Ky, W] = meshgrid(ky, kx, w);
K = sqrt(Kx.^2 + Ky.^2);

%% Windowing Function
window_xy = hann(Nx) * hann(Ny)';
window_t = hann(Nt);
window = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% Pre-allocate results
s_Date = NaT(length(file_list), 1);
s_Tp_surf = zeros(length(file_list), 1);
s_Tp_wave = zeros(length(file_list), 1);
s_Pdir_surf = zeros(length(file_list), 1);
s_Pdir_wave = zeros(length(file_list), 1);
s_Signal_surf = zeros(length(file_list), 1);
s_Signal_wave = zeros(length(file_list), 1);
s_Noise_surf = zeros(length(file_list), 1);
s_Noise_wave = zeros(length(file_list), 1);
s_Land_energy = zeros(length(file_list), 1);

%% Start parallel loop
tic
for i = 1:1%length(file_list)

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

    %% Windowing
    png_wave_no = png_wave;
    png_wave = png_wave .* window;

    %% FFT (with spectrum normalization)
    image_spectra_wave_no = abs(fftn(png_wave_no)).^2 / (Lx * Ly * Lt);
    image_spectra_wave_no = fftshift(image_spectra_wave_no);

    image_spectra_wave = abs(fftn(png_wave)).^2 / (Lx * Ly * Lt);
    image_spectra_wave = fftshift(image_spectra_wave);

    %% High-Pass Filter
    W_pass = 0.35;
    K_pass = 0.03;
    filter_HP = (abs(W) < W_pass) .* (K <= K_pass);
    filter_HP = ~ filter_HP;

    image_spectra_wave_HP = image_spectra_wave .* filter_HP;

    %% 2D Image Spectrum
    filter_2D = (W > 0);

    image_spectra_wave_2D = image_spectra_wave_HP .* filter_2D;
    image_spectra_wave_2D = 2 * image_spectra_wave_2D;

    %% Current Velocity
    E = image_spectra_wave_2D .* K.^(1.21);

    sigma = sqrt(g * K .* tanh(K * h));

    E_Kx_Kx = sum(E .* Kx.^2, 'all');
    E_Ky_Ky = sum(E .* Ky.^2, 'all');
    E_Kx_Ky = sum(E .* Kx .* Ky, 'all');
    E_W_Sigma_Kx_Kx = sum(E .* (W - sigma) .* Kx.^2, 'all');
    E_W_Sigma_Ky_Ky = sum(E .* (W - sigma) .* Ky.^2, 'all');

    wave_ux = (E_Ky_Ky * E_W_Sigma_Kx_Kx - E_Kx_Ky * E_W_Sigma_Ky_Ky) / (E_Kx_Kx * E_Ky_Ky - (E_Kx_Ky)^2);
    wave_uy = (E_Kx_Kx * E_W_Sigma_Ky_Ky - E_Kx_Ky * E_W_Sigma_Kx_Kx) / (E_Kx_Kx * E_Ky_Ky - (E_Kx_Ky)^2);

    wave_U = [wave_ux, wave_uy];
    disp(wave_U);

    %% Band-Pass Filter
    BPV =  8 * (2 * pi / Lt);

    filter_wave_BP = (abs(W - (sigma + Kx * wave_ux + Ky * wave_uy)) < BPV) | (abs(W - (-sigma + Kx * wave_ux + Ky * wave_uy)) < BPV);

    image_spectra_wave_BP = image_spectra_wave .* filter_wave_BP;    

    %% Modulation Transfer Function
    image_spectra_wave_final = image_spectra_wave_BP .* K.^(1.21);

    %% Signal to Noise Ratio
    signal_wave = sum(image_spectra_wave_final, 'all');

    noise_wave = sum(image_spectra_wave_no, 'all') - signal_wave;


    %% Save results
    s_Signal_wave(i) = signal_wave;
    s_Noise_wave(i) = noise_wave;
    s_Date(i) = date;

    %% Checker
    disp([length(file_list), i]);

end
toc


save("SNR01.mat", 's_Date', 's_Signal_wave', 's_Noise_wave');