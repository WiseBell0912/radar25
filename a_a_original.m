clear; close all; clc;

file_path = 'E:/png2019/10/';
flist = dir([file_path, '*.png']);

%% save result

for ii = 1:length(flist)
    pngname = [file_path, flist(ii).name]
    [SNR(ii), signal(ii), noise(ii), Tp(ii), Pdir(ii)] = cal1st(pngname);
end

save('result1st_surf.mat','SNR','signal','noise','Tp','Pdir')

function [SNR, signal, noise, Tp, Pdir] = cal1st(pngname)
%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(165);
surf_center = 830 + 630/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 630/2;

%% Radar Sub-Image sequence
date = pngname(end-16:end-4);
date = datetime(date, 'InputFormat', 'yyyyMMdd_HHmm');

png_long = imread(pngname);
png_long = png_long(1 : 512 * 1080 * 128);
png_long = reshape(png_long, 512, 1080, 128);
png_long = flip(png_long, 2);
png_long = flip(png_long, 3);

[png_surf, png_wave] = f_extract_zone_old(png_long, modifiy_theta, surf_theta1, surf_theta2, surf_center, wave_theta1, wave_theta2, wave_center);

img_cut = png_wave;

%% input parameter
Nx = size(img_cut,2);
Ny = size(img_cut,1);
Nt = size(img_cut,3);

dx = 3;
dy = 3;
dt = 1.43;

Lx = Nx * dx - dx;
Ly = Ny * dy - dy;
Lt = Nt * dt - dt;

g = 9.81;

h = 50;

%% fourier transform

wind2d = hann(Ny) * hann(Nx)';

windT = hann(Nt);
wind = repmat(wind2d, 1, 1, Nt) .* reshape(windT, 1, 1, Nt);

ft = fftn(wind.*img_cut);%rad_img & waves
ft_ = fftn(img_cut);% inverse fft에 사용
% ft = fftshift(ft);% ft는 fftshift 사용x / pt만 fftshift

w_ft = 0 : 2*pi/Lt : 2*pi*(Nt-1)/Lt; w_ft = w_ft - w_ft(round(end/2)+1);
kx_ft = 0 : 2*pi/Lx : 2*pi*(Nx-1)/Lx; kx_ft = kx_ft - kx_ft(round(end/2)+1);
ky_ft = 0 : 2*pi/Ly : 2*pi*(Ny-1)/Ly; ky_ft = ky_ft - ky_ft(round(end/2)+1);

[Kx, Ky,W] = meshgrid(kx_ft,ky_ft,w_ft);
Kx = single(Kx); Ky = single(Ky); W = single(W);
K = sqrt(Kx.^2 + Ky.^2);

%% High pass filter
idxK = K(:,:,1) > 0.03;            idxK = fftshift(idxK);
idxw = abs(w_ft) > 0.35;    idxw = fftshift(idxw);

for ii = 1:Nt
    ft(:,:,ii) = ft(:,:,ii).*idxK * idxw(ii);
    ft_(:,:,ii) = ft_(:,:,ii).*idxK * idxw(ii);
end

pt_hp = abs(ft).^2 / Nt^2 / Nx^2 / Ny^2;
pt_hp_ = abs(ft_).^2 / Nt^2 / Nx^2 / Ny^2;% for SNR(jj)
%% current estimation
%%%%%%%%%%%%% full frequency algorithm
Kx1d = reshape(Kx,1,[]);    Ky1d = reshape(Ky,1,[]);
K1d =  sqrt(Kx1d.^2 + Ky1d.^2);     sigma1d = sqrt(g*K1d .* tanh(K1d*h));

pt_hp = fftshift(pt_hp);

var = 0;%var = [www sigma kxx kyy kkk ppp*kkk^(-1.2)*abs(www)^(-0.6)];
for ii = Nt/2:Nt
    pt1d = reshape(pt_hp(:,:,ii),1,[]);
    idx = pt1d > 0;%0.00001*max(max(max(pt)));

    Npt = sum(idx);
    Nend = size(var,1);

    if Npt ~=0
        var(Nend+1:Nend+Npt,1) = w_ft(ii);
        var(Nend+1:Nend+Npt,2) = sigma1d(idx);
        var(Nend+1:Nend+Npt,3) = Kx1d(idx);
        var(Nend+1:Nend+Npt,4) = Ky1d(idx);
        var(Nend+1:Nend+Npt,5) = K1d(idx);

        MTF1d = sqrt(K1d(idx).^(-1.2) * abs(w_ft(ii))^(-0.6));

        var(Nend+1:Nend+Npt,6) = pt1d(idx) .* MTF1d;
    end

end

var = var(2:end,:);

%%%%%%%%%%%%% full frequency algorithm


idx = var(:,1)<0;
var(idx,2) = -var(idx,2);

a1 = sum(var(:,6).*var(:,3).^2);%kx^2
b1 = sum(var(:,6).*var(:,3).*var(:,4));%kx*ky
c1 = sum(var(:,6).*(var(:,1) - var(:,2)) .* var(:,3));%(w-sigma)kx
a2 = b1;%kx*ky
b2 = sum(var(:,6).*var(:,4).^2);%ky^2
c2 = sum(var(:,6).*(var(:,1) - var(:,2)) .* var(:,4));%(w-sigma)ky

if (a1*b2 - a2*b1) ~= 0
    ux_cal = (c1*b2 - c2*b1) / (a1*b2 - a2*b1);
    uy_cal = (c1*a2 - c2*a1) / (a1*b2 - a2*b1);
else
    ux_cal = a1/c1;
    uy_cal = a2/c2;
    disp('error or 1d wave')
    ux
    uy
end

disp(ux_cal);
disp(uy_cal);

clear var

%% band pass filter
wk1 = sqrt(g*K.*tanh(K*h)) + Kx*ux_cal + Ky*uy_cal;
wk2 = -sqrt(g*K.*tanh(K*h)) + Kx*ux_cal + Ky*uy_cal;

bpv = 2 * (Nt/32) * 2*pi/Lt; %Band pass filtering 하는 정도 / a * (dw)
idx = ~((W > wk1 - bpv & W < wk1 + bpv) | (W > wk2 - bpv & W < wk2 + bpv));
idx = fftshift(idx);

ft(idx) = 0;
ft_(idx) = 0;
pt_bp = abs(ft).^2 / Nt^2 / Nx^2 / Ny^2;
pt_bp_ = abs(ft_).^2 / Nt^2 / Nx^2 / Ny^2;%windowing 미적용 for SNR(jj)

%% spectrum analysis
% pt = fftshift(pt);
pt_bp_ = fftshift(pt_bp_);

Sk2d_bp = (sum(pt_bp_(:,:,round(end/2):end),3));

pw1d = 2*squeeze(sum(sum(pt_bp_,1),2));
pw1d(1:end/2) = 0;

%% calculate SNR
signal = sum(sum(sum(pt_bp_)));
noise = sum(sum(sum(pt_hp_))) - signal;
SNR = sqrt(signal / noise);

%% calculate Tp
[vv, idx] = max(pw1d);

Tp = 2*pi/w_ft(idx);

%% calculate Pdir
idx = Sk2d_bp == max(max(Sk2d_bp));

Pdir = mod((90 - rad2deg(atan2(Ky(idx),Kx(idx))) ) - rad2deg(beta) - mdeg,360);
end

