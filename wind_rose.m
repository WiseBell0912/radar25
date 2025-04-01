clear; close all; clc;

%% 불러오기
% ADCP
load("ADCP.mat");
a_Date = a_Date02;
a_Hs = a_Hs02;
a_Pdir = a_Pdir02;
a_Tp = a_Tp02;

clear a_Date01 a_Hs01 a_Pdir01 a_Tp01 a_Date02 a_Hs02 a_Pdir02 a_Tp02

% Bouy
load("new_Bouy_no_NAN_ZERO.mat");

% Radar wave
load("./0327/snr_y1910_wave_0327.mat");
r_Date = Date;
r_wave_Pdir = Pdir;
r_wave_SNR = SNR;
r_wave_Tp = Tp;
r_wave_Ux = Ux;
r_wave_Uy = Uy;

load("./0327/snr_y1911_wave_0327.mat");
r_Date = [r_Date ; Date];
r_wave_Pdir = [r_wave_Pdir ; Pdir];
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_Tp = [r_wave_Tp ; Tp];
r_wave_Ux = [r_wave_Ux ; Ux];
r_wave_Uy = [r_wave_Uy ; Uy];

load("./0327/snr_y1912_wave_0327.mat");
r_Date = [r_Date ; Date];
r_wave_Pdir = [r_wave_Pdir ; Pdir];
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_Tp = [r_wave_Tp ; Tp];
r_wave_Ux = [r_wave_Ux ; Ux];
r_wave_Uy = [r_wave_Uy ; Uy];

r_wave_SNR = sqrt(r_wave_SNR);

% Radar surf
load("./0327/snr_y1910_surf_0327.mat");
r_Date = Date;
r_surf_Pdir = Pdir;
r_surf_SNR = SNR;
r_surf_Tp = Tp;
r_surf_Ux = Ux;
r_surf_Uy = Uy;

load("./0327/snr_y1911_surf_0327.mat");
r_Date = [r_Date ; Date];
r_surf_Pdir = [r_surf_Pdir ; Pdir];
r_surf_SNR = [r_surf_SNR ; SNR];
r_surf_Tp = [r_surf_Tp ; Tp];
r_surf_Ux = [r_surf_Ux ; Ux];
r_surf_Uy = [r_surf_Uy ; Uy];

load("./0327/snr_y1912_surf_0327.mat");
r_Date = [r_Date ; Date];
r_surf_Pdir = [r_surf_Pdir ; Pdir];
r_surf_SNR = [r_surf_SNR ; SNR];
r_surf_Tp = [r_surf_Tp ; Tp];
r_surf_Ux = [r_surf_Ux ; Ux];
r_surf_Uy = [r_surf_Uy ; Uy];

r_surf_SNR = sqrt(r_surf_SNR);

clear Date SNR Ux Uy r_wave_Ux r_wave_Uy r_surf_Ux r_surf_Uy

%% 처리 
mask_r2a = ismember(r_Date, a_Date);

r2a_Date = r_Date(mask_r2a);
r2a_wave_Pdir = r_wave_Pdir(mask_r2a);
r2a_wave_Tp = r_wave_Tp(mask_r2a);
r2a_surf_Pdir = r_surf_Pdir(mask_r2a);
r2a_surf_Tp = r_surf_Tp(mask_r2a);

mask_a2r = ismember(a_Date, r_Date);

a2r_Date = a_Date(mask_a2r);
a2r_Pdir = a_Pdir(mask_a2r);
a2r_Tp = a_Tp(mask_a2r);

%% 그래프 Pdir
figure(1)
tiledlayout(1,3);

nexttile;
polarhistogram(deg2rad(r2a_wave_Pdir), 72);
title('Wave Pdir');

nexttile;
polarhistogram(deg2rad(a2r_Pdir), 72);
title('ADCP Pdir');

nexttile;
polarhistogram(deg2rad(r2a_surf_Pdir), 72);
title('Surf Pdir');

%% 그래프 Pdir
figure(2)
tiledlayout(2,1);

nexttile;
hold on;
plot(r2a_Date, r2a_wave_Pdir-180.*(r2a_wave_Pdir > 180));
plot(r2a_Date, r2a_surf_Pdir-180.*(r2a_surf_Pdir > 180));
plot(a2r_Date,a2r_Pdir-180.*(a2r_Pdir > 180));
hold off;
title('Pdir');
legend('Wave', 'Surf', 'ADCP');
xlim([datetime(2019, 10, 1), datetime(2019, 11, 1)])

%% 그래프 Tp
nexttile;

hold on;
plot(r2a_Date, r2a_wave_Tp);
plot(r2a_Date, r2a_surf_Tp);
plot(a2r_Date, a2r_Tp);
hold off;
title('Tp');
legend('Wave', 'Surf', 'ADCP');
xlim([datetime(2019, 10, 1), datetime(2019, 11, 1)]);