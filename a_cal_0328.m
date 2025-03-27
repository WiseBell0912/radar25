clear; close all; clc;

%% Search
file_path = 'E:/png2019/11/';
file_list = dir([file_path, '*.png']);

%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(160);
surf_center = 870 + 600/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 600/2;

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
h = 50;

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
sc        = Nt * Nx * Ny / sum(window(:).^2);

%% 결과 저장용
nFile = length(file_list);
Date  = NaT(nFile,1);
M0    = zeros(nFile,1);

%%
tic;

for i = 1 : nFile
    %% 파일 읽기
    png_path = fullfile(file_list(i).folder, file_list(i).name);
    dateStr  = file_list(i).name(5:end-4);
    dateVal  = datetime(dateStr, 'InputFormat', 'yyyyMMdd_HHmm');

    png_long = imread(png_path);
    png_long = png_long(1 : 512*1080*128);  % 파일 크기에 맞춰 잘라내는 로직(원 코드 동일)
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    [png_surf, png_wave] = f_extract_zone( ...
        png_long, modifiy_theta, ...
        surf_theta1, surf_theta2, surf_center, ...
        wave_theta1, wave_theta2, wave_center);

    img = single(png_wave);
    png_long = [];
    png_surf = [];
    png_wave = [];

    %% FFT
    img_windowed          = img .* window;
    image_spectrum        = fftn(img_windowed);

    %% M0
    m0Val = sum(image_spectrum, 'all');

    %% 결과 저장
    M0(i)  = m0Val;
    Date(i) = dateVal;

    i

end

toc;

save("snr_y1911_wave_0328.mat", 'Date', 'M0');
