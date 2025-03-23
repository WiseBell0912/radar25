clear; close all; clc;

%% Search
file_path = '/Volumes/임현종/png2019/10/';
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
for i = 1 : length(file_list)

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
    % need
    image = png_wave;
    sc = Nx * Ny * Nt / sum(window.^2, 'all');

    % windowing
    image_y_window = image .* window;
    image_n_window = image;

    % fft
    image_y_window_spectrum = fftshift(fftn(image_y_window));
    image_n_window_spectrum = fftshift(fftn(image_n_window));

    % spectrum
    image_y_window_spectrum = abs(image_y_window_spectrum).^2 / Nx^2 / Ny^2 / Nt^2;
    image_n_window_spectrum = abs(image_n_window_spectrum).^2 / Nx^2 / Ny^2 / Nt^2;
    pt = sc*abs(image_y_window_spectrum).^2 / Nx^2 / Ny^2 / Nt^2;

    % hpf
    fK = K > 0.03;
    fW = W > 0.35;
    image_y_window_spectrum_hpf = image_y_window_spectrum .* fK .* fW;
    image_n_window_spectrum_hpf = image_n_window_spectrum .* fK .* fW;

    % current velocity
    fCV = (image_y_window_spectrum_hpf > 0.00001 * max(pt, [], 'all') & W > 0);
    sigma = sqrt(g .* K .* tanh(K .* h));
    weight = image_y_window_spectrum_hpf;% .* sqrt(K.^(-1.2) .* abs(wave_center).^(0.6));

    a1 = weight .* Kx .* Kx .* fCV;
    b1 = weight .* Kx .* Ky .* fCV;
    c1 = weight .* (W - sigma) .* Kx .* fCV;
    a2 = weight .* Kx .* Ky .* fCV;
    b2 = weight .* Ky .* Ky .* fCV;
    c2 = weight .* (W - sigma) .* Ky .* fCV;

    a1 = sum(a1, 'all');
    b1 = sum(b1, 'all');
    c1 = sum(c1, 'all');
    a2 = sum(a2, 'all');
    b2 = sum(b2, 'all');
    c2 = sum(c2, 'all');

    if (a1*b2 - a2*b1) ~= 0
        ux_cal = (c1*b2 - c2*b1) / (a1*b2 - a2*b1)
        uy_cal = (c1*a2 - c2*a1) / (a1*b2 - a2*b1)
    else
        ux_cal = a1/c1
        uy_cal = a2/c2
        disp('error or 1d wave')
        ux
        uy
    end


    %% Checker
    disp([length(file_list), i]);

end
toc