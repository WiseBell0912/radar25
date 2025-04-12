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
load("Bouy.mat");
b_Pdir = b_WaveDirection;
b_Tp = b_MaximumWavePeriod;

% Radar wave
load("snr_y1910_0411.mat");
r_WaveE = WaveE;
r_wave_Uy = wave_Uy;
r_wave_Ux = wave_Ux;
r_wave_Tp = wave_Tp;
r_wave_SNR = wave_SNR;
r_wave_Pdir = wave_Pdir;
r_SurfE = SurfE;
r_surf_Uy = surf_Uy;
r_surf_Ux = surf_Ux;
r_surf_Tp = surf_Tp;
r_surf_SNR = surf_SNR;
r_surf_Pdir = surf_Pdir;
r_LandE = LandE;
r_Date = Date;

load("snr_y1911_0411.mat");
r_WaveE = [r_WaveE; WaveE];
r_wave_Uy = [r_wave_Uy; wave_Uy];
r_wave_Ux = [r_wave_Ux; wave_Ux];
r_wave_Tp = [r_wave_Tp; wave_Tp];
r_wave_SNR = [r_wave_SNR; wave_SNR];
r_wave_Pdir = [r_wave_Pdir; wave_Pdir];
r_SurfE = [r_SurfE; SurfE];
r_surf_Uy = [r_surf_Uy; surf_Uy];
r_surf_Ux = [r_surf_Ux; surf_Ux];
r_surf_Tp = [r_surf_Tp; surf_Tp];
r_surf_SNR = [r_surf_SNR; surf_SNR];
r_surf_Pdir = [r_surf_Pdir; surf_Pdir];
r_LandE = [r_LandE; LandE];
r_Date = [r_Date; Date];

load("snr_y1912_0411.mat");
r_WaveE = [r_WaveE; WaveE];
r_wave_Uy = [r_wave_Uy; wave_Uy];
r_wave_Ux = [r_wave_Ux; wave_Ux];
r_wave_Tp = [r_wave_Tp; wave_Tp];
r_wave_SNR = [r_wave_SNR; wave_SNR];
r_wave_Pdir = [r_wave_Pdir; wave_Pdir];
r_SurfE = [r_SurfE; SurfE];
r_surf_Uy = [r_surf_Uy; surf_Uy];
r_surf_Ux = [r_surf_Ux; surf_Ux];
r_surf_Tp = [r_surf_Tp; surf_Tp];
r_surf_SNR = [r_surf_SNR; surf_SNR];
r_surf_Pdir = [r_surf_Pdir; surf_Pdir];
r_LandE = [r_LandE; LandE];
r_Date = [r_Date; Date];

r_wave_SNR = sqrt(r_wave_SNR);
r_surf_SNR = sqrt(r_surf_SNR);

clear Date SNR Ux Uy r_wave_Ux r_wave_Uy r_surf_Ux r_surf_Uy

%% 처리 
mask1 = ismember(r_Date, a_Date);
mask2 = ismember(r_Date, b_Date);
mask = mask1 & mask2;

rr_Date = r_Date(mask);
rr_wave_Pdir = r_wave_Pdir(mask);
rr_wave_Tp = r_wave_Tp(mask);
rr_surf_Pdir = r_surf_Pdir(mask);
rr_surf_Tp = r_surf_Tp(mask);

mask1 = ismember(a_Date, b_Date);
mask2 = ismember(a_Date, r_Date);
mask = mask1 & mask2;

aa_Date = a_Date(mask);
aa_Pdir = a_Pdir(mask);
aa_Tp = a_Tp(mask);

mask1 = ismember(b_Date, a_Date);
mask2 = ismember(b_Date, r_Date);
mask = mask1 & mask2;

bb_Date = b_Date(mask);
bb_Pdir = b_Pdir(mask);
bb_Tp = b_Tp(mask);
%% 그래프 Pdir
figure(1)
tiledlayout(1,4);

nexttile;
polarhistogram(deg2rad(rr_wave_Pdir), 72);
title('Wave Pdir');

nexttile;
polarhistogram(deg2rad(rr_surf_Pdir), 72);
title('Surf Pdir');

nexttile;
polarhistogram(deg2rad(aa_Pdir), 72);
title('ADCP Pdir');

nexttile;
polarhistogram(deg2rad(bb_Pdir), 72);
title('Bouy Pdir');

%% 그래프 Pdir
figure(2)
tiledlayout(2,1);

nexttile;
hold on;
% plot(rr_Date, rr_wave_Pdir-180.*(rr_wave_Pdir > 180));
% plot(rr_Date, rr_surf_Pdir-180.*(rr_surf_Pdir > 180));
% plot(aa_Date,aa_Pdir-180.*(aa_Pdir > 180));
% plot(aa_Date,bb_Pdir-180.*(bb_Pdir > 180));
plot(rr_Date, rr_wave_Pdir);
plot(rr_Date, rr_surf_Pdir);
plot(aa_Date,aa_Pdir);
plot(aa_Date,bb_Pdir);
hold off;
title('Pdir');
legend('Wave', 'Surf', 'ADCP', 'Bouy');
xlim([datetime(2019, 10, 1), datetime(2019, 11, 1)])

%% 그래프 Tp
nexttile;

hold on;
plot(rr_Date, rr_wave_Tp);
plot(rr_Date, rr_surf_Tp);
plot(aa_Date, aa_Tp);
plot(bb_Date, bb_Tp);
hold off;
title('Tp');
legend('Wave', 'Surf', 'ADCP', 'Bouy');
xlim([datetime(2019, 10, 1), datetime(2019, 11, 1)]);