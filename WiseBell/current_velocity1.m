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

%% Frequency
kx = -pi/dx : 2*pi/Lx : pi/dx;
ky = -pi/dy : 2*pi/Ly : pi/dy;
w = -pi/dt : 2*pi/Lt : pi/dt;
w(65) = [];

[Kx2, Ky2] = meshgrid(kx, ky);
K2 = sqrt(Kx2.^2 + Ky2.^2);

[Kx3, Ky3, W3] = meshgrid(kx, ky, w);
K3 = sqrt(Kx3.^2 + Ky3.^2);

%% Dispersion shell
sigma = sqrt(g.*K2.*tanh(K2.*h));
figure(1);
surf(Kx2, Ky2, sigma, 'EdgeAlpha', 0, 'FaceAlpha', 0.9);
title("Dispersion Shell without Surface Current");
xlabel("kx"); ylabel("ky"); zlabel("w");
xlim([-2, 2]); ylim([-2, 2]);

dispersion = sqrt(g.*K2.*tanh(K2.*h)) + Kx2.*0 + Ky2.*2;
figure(2);
surf(Kx2, Ky2, dispersion, 'EdgeAlpha', 0, 'FaceAlpha', 0.9);
title("Dispersion Shell with Surface Current");
xlabel("kx"); ylabel("ky"); zlabel("w");
xlim([-2, 2]); ylim([-2, 2]);

%% Windowing Function
window_xy = hann(Ny) * hann(Nx)';
window_t = hann(Nt);
window = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% BP
K_BP = ones(211, 121, 130);

%% Save
Date = NaT(length(file_list), 1);
Ux = zeros(length(file_list), 1);
Uy = zeros(length(file_list), 1);

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

    %% Windowing
    png_wave = png_wave .* window;

    %% FFT (with spectrum normalization)
    image_spectra_wave = abs(fftn(png_wave)).^2 / Nx / Ny / Nt;
    image_spectra_wave = fftshift(image_spectra_wave);

    %% Surface Current
    E = 2 .* sum(image_spectra_wave, 3) .* K2.^(1.21).* (W3 >0);

    A11 = sum(E .* Kx3.^2, "all");
    A12 = sum(E .* Kx3 .* Ky3, 'all');
    A21 = A12;
    A22 = sum(E .* Ky3.^2, 'all');

    B1 = sum(E .* (W3 - sqrt(g.*K3.*tanh(K3.*h)).*Kx3), 'all');
    B2 = sum(E .* (W3 - sqrt(g.*K3.*tanh(K3.*h)).*Ky3), 'all');

    A = [A11, A12; A21, A22];
    B = [B1; B2];

    U = A \ B;

    Ux_estimated = double(U(1));
    Uy_estimated = double(U(2));

    disp(['Estimated Ux: ', num2str(Ux_estimated), ' m/s']);
    disp(['Estimated Uy: ', num2str(Uy_estimated), ' m/s']);

    %% Save
    Date(i) = date;
    Ux(i) = Ux_estimated;
    Uy(i) = Uy_estimated;

    %% Checker
    disp([length(file_list), i]);

end
toc