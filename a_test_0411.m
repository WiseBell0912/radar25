clear; clc;

%% Search
file_path = 'E:/png2019/10/';
%file_path = 'E:/png2019/10/';
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

%% 하이패스 필터(HPF)용 마스크 미리 계산(반복문 밖)
hpW = (abs(W) > 2*pi*0.03);
hpMask = hpW;

%% 결과 저장용
nFile = length(file_list);

r_Date  = NaT(nFile,1);
r_LandE = zeros(nFile,1);
r_SurfE = zeros(nFile,1);
r_WaveE = zeros(nFile,1);

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

for i = 785 % 반절
%for i = 685 % 전체
%for i = 2564 % 전체
    %% 파일 읽기 및 Zone 추출
    png_path = fullfile(file_list(i).folder, file_list(i).name);
    dateStr  = file_list(i).name(5:end-4);
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

    img_surf = rot90(flip(flip(img_surf), 2));

    img_spectrum = fftn(img_surf .* window);
    img_spectrum = abs(img_spectrum).^2 / Nx^2 / Ny^2 / Nt^2;
    img_spectrum = fftshift(img_spectrum) .* hpW;

    % Current
    sigma_surf = sqrt(g .* K .* tanh(K .* 15));

    MTF_surf = (K.^(-1.21)) .* (W > 0) .* (img_spectrum > 1e-7 * max(img_spectrum(:)));
    MTF_surf(~isfinite(MTF_surf)) = 0;

    a1 = sum( img_spectrum .* MTF_surf .* Kx.^2 , 'all');
    a2 = sum( img_spectrum .* MTF_surf .* Ky.^2 , 'all');
    a3 = sum( img_spectrum .* MTF_surf .* Kx .* Ky , 'all');
    a4 = sum( img_spectrum .* MTF_surf .* (W - sigma_surf) .* Kx , 'all');
    a5 = sum( img_spectrum .* MTF_surf .* (W - sigma_surf) .* Ky , 'all');

    ux_surf = ( a2 * a4 - a3 * a5 ) / ( a1 * a2 - a3^2 )
    uy_surf = ( a1 * a5 - a3 * a4 ) / ( a1 * a2 - a3^2 )

    % View
    for j = 1 : 128
        figure(1);
        tiledlayout(1, 2);

        nexttile;
        hold on;
        surf(img_surf(101:201, :, j), 'EdgeAlpha', 0);
        quiver3(101, 101, max(img_surf(:)) + 10, ux_surf, uy_surf, 0, 10, 'r', 'LineWidth', 2);
        hold off;
        view(0, 90);

        nexttile;
        hold on;
        surf(sum(img_spectrum, 3), 'EdgeAlpha', 0);
        quiver3(101, 101, max(img_spectrum(:)) + 10, ux_surf, uy_surf, 0, 10, 'r', 'LineWidth', 2);
        hold off;
        view(0, 90);

    end

end