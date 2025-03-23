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

rdata = imread(pngname);

%%

%%%%%%%%%%%%%%
sdrng = 800; %[m]
dr = 3; %[m]
nr = 512;
mdeg = 76;% modification factor for north correction

rr = sdrng:dr:sdrng+dr*(nr-1);

nimg = 128;
nbearing = 1080;

theta = linspace(0,2*pi,nbearing);

[Th,R] = meshgrid(theta,rr);

X = R.*cos( deg2rad(90) - Th - deg2rad(mdeg));
Y = R.*sin( deg2rad(90) - Th - deg2rad(mdeg));
%%%%%%%%%%%%%%
rdata = rdata(1:nr*nbearing*nimg);
rdata = reshape(rdata,nr,nbearing,nimg);
rdata = flip(rdata,3);
%% pol2car

%%%%%%%%%%%%% cut info 
alpha_deg = 329;% box position rotation degree
% beta_deg = 329;% box rotation degree
mdeg = 76;% modification factor for north correction
Ly = 630;%radial direction length / Ly
Lx = 360;%perpendicular direction length / Lx
rcenter = 800+630/2;% radar center to box center range

alpha = deg2rad(90 - alpha_deg);
beta = alpha;

%%%%%%%%%%%%% cut info

% divide radar image
%%%%%%%%%%%%%%%%%%%%%%%%% cut process
%%%%%%%%%% cut 
x = -Lx/2:dr:Lx/2;
y = -Ly/2:dr:Ly/2;
[Xc,Yc] = meshgrid(x,y);

Xtemp = Yc*cos(beta) - Xc*sin(beta);
Ytemp = Yc*sin(beta) + Xc*cos(beta);

Xc = Xtemp + rcenter*cos(alpha);
Yc = Ytemp + rcenter*sin(alpha);

img_cut = single(zeros(size(Xc,1),size(Xc,2),nimg));

% cut range for plot
cutx = Lx*[-0.5 0.5 0.5 -0.5 -0.5];
cuty = Ly*[0.5 0.5 -0.5 -0.5 0.5];

cutxx = cuty*cos(beta) - cutx*sin(beta);
cutyy = cuty*sin(beta) + cutx*cos(beta);

cutX = cutxx + rcenter*cos(alpha);
cutY = cutyy + rcenter*sin(alpha);

%%%%%%%%%% cut 

%%% idx mapping
temp = Xc + Yc*1i;
temp_r = abs(temp);
temp_t = angle(temp);
temp_t = mod( pi/2 - temp_t - deg2rad(mdeg),2*pi);

for ii = 1:nimg
    img_temp = rdata(:,:,ii);
    
    idx_Tc = floor(temp_t/abs(theta(2)-theta(1))) + 1;
    idx_Rc = floor((temp_r - sdrng)/dr) + 1;
    idx = (idx_Rc + nr * (idx_Tc - 1));
    temp = img_temp(idx);
    img_cut(:,:,ii) = temp;
    
end
%%%%%%%%%% cut 


%% input parameter
Nx = size(img_cut,2);
Ny = size(img_cut,1);
Nt = size(img_cut,3);

dt = 1.43;

Lt = Nt*dt;

h = 30;%water depth
g = 9.81;

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

