clear; close all; clc;

%% Search
file_path = 'D:/png2019/10/';
file_list = dir([file_path, '*.png']);

%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(165);
surf_center = 830 + 630/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 630/2;

%% Save
Date = NaT(length(file_list), 1);
Wave_E = zeros(length(file_list), 1);
Surf_E = zeros(length(file_list), 1);

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

    Date(i) = date;
    Wave_E(i) = sum(png_wave, 'all');
    Surf_E(i) = sum(png_surf, 'all');


    %% Checker
    disp([length(file_list), i]);

end
toc