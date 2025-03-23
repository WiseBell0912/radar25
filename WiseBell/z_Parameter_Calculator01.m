%% 해상 파라미터를 추출하는 스크립트
clear; close all; clc;

%% 레이더 이미지 확인
path = 'C:/Users/Hyeonjong Im/Documents/GitHub/Radar_2nd/PNG/';
file_list = dir([path, '*.png']);

%% 레이더 이미지 만큼 반복하며 해상 파라미터 추출

for i = 1 : 
    png_location = [path, file_list(i).name];
    [Date(i), surf_SNR(i), surf_Signal(i), surf_Noise(i), surf_Tp(i), wave_SNR(i), wave_Signal(i), wave_Noise(i), wave_Tp(i), wave_Pdir(i)] = cal_wave_para(png_location, j);

    disp(i); disp('/'); disp(length(file_list));
end

%% 추출한 해상 파라미터 저장
%save();

%% 해상 파라미터 추출 함수
function [Date, surf_SNR, surf_Signal, surf_Noise, surf_Tp, wave_SNR, wave_Signal, wave_Noise, wave_Tp, wave_Pdir] = cal_wave_para(png_location, h)

% 레이더 이미지 불러오기
%png_location = ['/Users/imhyeonjong/Documents/POL/', file_list(1).name];
png = imread(png_location);

% 날짜
Date = datetime(png_location(end-16:end-4), 'Inputformat', 'yyyyMMdd_HHmm');

% 레이더 이미지 기본 정보
png_count = 128;

r_count = 512;
r_in = 800;
r_out = 2333;
r_delta = 3;

th_count = 1080;
th_modify = 62;

% pngLong 을 위한 그리드 생성
r = r_in : r_delta : r_out;

th = linspace(0, 2*pi, th_count);

[Th, R] = meshgrid(th, r);

X = R .* cos(Th + deg2rad(th_modify));
Y = R .* sin(Th + deg2rad(th_modify));

% 타임 시퀀스에 맞도록 3차원으로 변경
png = png(1 : r_count*th_count*png_count);
png = reshape(png, r_count, th_count, png_count);
png = flip(png, 3);

% 실제와 같게 이미지 좌우 반전
png = flip(png, 2);

%% Surf Box
% Surf Box 기본 정보
box_deg = (90 - th_modify + 14) - 121;
box_rad = deg2rad(box_deg);
box_x = 630;
box_y = 360;
box_center = 900 + box_x/2;

% Surf Box 그리드 생성
x = -box_x/2 : r_delta : box_x/2;
y = -box_y/2 : r_delta : box_y/2;
[Xc,Yc] = meshgrid(x,y);

Xtemp = Xc*cos(box_rad) - Yc*sin(box_rad);
Ytemp = Xc*sin(box_rad) + Yc*cos(box_rad);

Xc = Xtemp + box_center*cos(box_rad);
Yc = Ytemp + box_center*sin(box_rad);

surf_GridX = Xc;
surf_GridY = Yc;

% Surf Box 이미지 추출 및 저장
img_cut_surf = single( zeros(size(Xc,1), size(Xc,2), png_count) );

temp = Xc + Yc * 1i;
temp_r = abs(temp);
temp_th = angle(temp);
temp_th = mod(deg2rad(90 - th_modify) - temp_th, 2*pi);

for ii = 1 : png_count
    img_temp = png(:,:,ii);
    
    idx_Tc = floor(temp_th/abs(th(2)-th(1))) + 1;
    idx_Rc = floor((temp_r - r_in)/r_delta) + 1;
    idx = (idx_Rc + r_count * (idx_Tc - 1));
    temp = img_temp(idx);
    img_cut_surf(:,:,ii) = temp;
end

surf_Png = img_cut_surf;

%% Wave Box
% Wave Box 기본 정보
box_deg = (90 - th_modify + 14) - 84;
box_rad = deg2rad(box_deg);
box_x = 630;
box_y = 360;
box_center = 1150 + box_x/2;

% Wave Box 그리드 생성
x = -box_x/2 : r_delta : box_x/2;
y = -box_y/2 : r_delta : box_y/2;
[Xc,Yc] = meshgrid(x,y);

Xtemp = Xc*cos(box_rad) - Yc*sin(box_rad);
Ytemp = Xc*sin(box_rad) + Yc*cos(box_rad);

Xc = Xtemp + box_center*cos(box_rad);
Yc = Ytemp + box_center*sin(box_rad);

wave_GridX = Xc;
wave_GridY = Yc;

% Wave Box 이미지 추출 및 저장
img_cut_wave = single( zeros(size(Xc,1), size(Xc,2), png_count) );

temp = Xc + Yc * 1i;
temp_r = abs(temp);
temp_th = angle(temp);
temp_th = mod(deg2rad(90 - th_modify) - temp_th, 2*pi);

for ii = 1 : png_count
    img_temp = png(:,:,ii);
    
    idx_Tc = floor(temp_th/abs(th(2)-th(1))) + 1;
    idx_Rc = floor((temp_r - r_in)/r_delta) + 1;
    idx = (idx_Rc + r_count * (idx_Tc - 1));
    temp = img_temp(idx);
    img_cut_wave(:,:,ii) = temp;
    
end

wave_Png = img_cut_wave;

%% 기본 물리적 정보
Nx = size(surf_Png, 2);
Ny = size(surf_Png, 1);
Nt = size(surf_Png, 3);

dx = 3;
dy = 3;
dt = 1.43;

Lx = 630;
Ly = 360;
Lt = dt * Nt;

h = 30;
g = 9.81;

%% Image normalization
surf_Png = surf_Png - mean(surf_Png);
wave_Png = wave_Png - mean(wave_Png);

%% Hanning window
wind_xy = hann(Ny) * hann(Nx)';
wind_t = hann(Nt);
wind = repmat(wind_xy, 1, 1, Nt) .* reshape(wind_t, 1, 1, Nt);

surf_Png = surf_Png .* wind;
wave_Png = wave_Png .* wind;

clearvars wind_xy wind_t;

%% Fourier transform
surf_power = abs(fftn(surf_Png)).^2;
wave_power = abs(fftn(wave_Png)).^2;

surf_power = fftshift(surf_power);
wave_power = fftshift(wave_power);

%% Power spectrum normalization
surf_power = surf_power ./ Nx ./ Ny ./ Nt;
wave_power = wave_power ./ Nx ./ Ny ./ Nt;

%% Grid setting
w = 0 : 2*pi/Lt : 2*pi*(Nt-1)/Lt;
kx = 0 : 2*pi/Lx : 2*pi*(Nx-1)/Lx;
ky = 0 : 2*pi/Ly : 2*pi*(Ny-1)/Ly;

w = w - w(round(end/2));
kx = kx - kx(round(end/2));
ky = ky - ky(round(end/2));

[Kx, Ky, W] = meshgrid(kx, ky, w);

Kx = single(Kx);
Ky = single(Ky);
W = single(W);

K = sqrt(Kx.^2 + Ky.^2);

%% High pass filter
W_pass = 2 * pi * 0.03;
filter_W = W_pass < abs(W);

K_pass = 0.011191;
filter_K = K_pass <= abs(K);

filter_HP = filter_W .* filter_K;

surf_power = surf_power .* filter_HP;
wave_power = wave_power .* filter_HP;

surf_power_HP = surf_power;
wave_power_HP = wave_power;

%% Current velocity estimation
sigma = sqrt(g .* K .* tanh(K .* h));

a1 = sum( surf_power .* Kx.^2 , 'all');
a2 = sum( surf_power .* Ky.^2 , 'all');
a3 = sum( surf_power .* Kx .* Ky , 'all');
a4 = sum( surf_power .* (W - sigma) .* Kx.^2 , 'all');
a5 = sum( surf_power .* (W - sigma) .* Ky.^2 , 'all');

b1 = sum( wave_power .* Kx.^2 , 'all');
b2 = sum( wave_power .* Ky.^2 , 'all');
b3 = sum( wave_power .* Kx .* Ky , 'all');
b4 = sum( wave_power .* (W - sigma) .* Kx.^2 , 'all');
b5 = sum( wave_power .* (W - sigma) .* Ky.^2 , 'all');

surf_ux = ( a2 * a4 - a3 * a5 ) / ( a1 * a2 - a3^2 );
surf_uy = ( a1 * a5 - a3 * a4 ) / ( a1 * a2 - a3^2 );

wave_ux = ( b2 * b4 - b3 * b5 ) / ( b1 * b2 - b3^2 );
wave_uy = ( b1 * b5 - b3 * b4 ) / ( b1 * b2 - b3^2 );

%% Band pass filter
disp_pass_surf_pos = sqrt( g .* K .* tanh( K .* h ) ) + (Kx .* surf_ux) + (Ky .* surf_uy);
disp_pass_surf_neg = -sqrt( g .* K .* tanh( K .* h ) ) + (Kx .* surf_ux) + (Ky .* surf_uy);

disp_pass_wave_pos = sqrt( g .* K .* tanh( K .* h ) ) + (Kx .* wave_ux) + (Ky .* wave_uy);
disp_pass_wave_neg = -sqrt( g .* K .* tanh( K .* h ) ) + (Kx .* wave_ux) + (Ky .* wave_uy);

disp_boundary = 8 * (2*pi/Lt);

filter_BP_surf_pos = (disp_pass_surf_pos - disp_boundary <= W) + (W <= disp_pass_surf_pos + disp_boundary);
filter_BP_surf_pos = filter_BP_surf_pos == 2;
filter_BP_surf_neg = (disp_pass_surf_neg - disp_boundary <= W) + (W <= disp_pass_surf_neg + disp_boundary);
filter_BP_surf_neg = filter_BP_surf_neg == 2;

filter_BP_wave_pos = (disp_pass_wave_pos - disp_boundary <= W) + (W <= disp_pass_wave_pos + disp_boundary);
filter_BP_wave_pos = filter_BP_wave_pos == 2;
filter_BP_wave_neg = (disp_pass_wave_neg - disp_boundary <= W) + (W <= disp_pass_wave_neg + disp_boundary);
filter_BP_wave_neg = filter_BP_wave_neg == 2;

filter_BP_surf = filter_BP_surf_pos + filter_BP_surf_neg;
filter_BP_wave = filter_BP_wave_pos + filter_BP_wave_neg;

surf_power = surf_power .* filter_BP_surf;
wave_power = wave_power .* filter_BP_wave;

%% SNR
surf_Signal = sum(surf_power, 'all');
wave_Signal = sum(wave_power, 'all');

surf_Noise = sum(surf_power_HP, 'all') - surf_Signal;
wave_Noise = sum(wave_power_HP, 'all') - wave_Signal;

surf_SNR = sqrt(surf_Signal / surf_Noise);
wave_SNR = sqrt(wave_Signal / wave_Noise);

%% MTF
surf_power = surf_power .* K.^(1.21);
wave_power = wave_power .* K.^(1.21);

%% Peak period
surf_frequency_spectrum = sum(sum(surf_power, 1), 2);
wave_frequency_spectrum = sum(sum(wave_power, 1), 2);

surf_frequency_spectrum = squeeze(surf_frequency_spectrum);
wave_frequency_spectrum = squeeze(wave_frequency_spectrum);

[dummy, idx_surf] = max(surf_frequency_spectrum);
[dummy, idx_wave] = max(wave_frequency_spectrum);

surf_Tp = 2*pi / abs(w(idx_surf));
wave_Tp = 2*pi / abs(w(idx_wave));

%% Wave direction
wave_directional_spectrum = sum(wave_power(end/2:end), 3);

[dummy, idx_wave] = max(wave_directional_spectrum, [], 'all');

wave_Pdir = mod( rad2deg(atan2(Ky(idx_wave), Kx(idx_wave))) - box_deg , 360);

end