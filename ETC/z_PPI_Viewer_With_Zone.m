clear; close all; clc;

% Year
list1 = {'2018', '2019', '2020', '2021', '2022', '2023'};
[index1, dummy1] = listdlg('PromptString', {'Select a Year.'}, 'SelectionMode', 'single', 'ListString', list1);
if dummy1 == 0
    return;
elseif strcmp(list1{index1}, '2018') == 1
    cd png2018
elseif strcmp(list1{index1}, '2019') == 1
    cd png2019
elseif strcmp(list1{index1}, '2020') == 1
    cd png2020
elseif strcmp(list1{index1}, '2021') == 1
    cd png2021
elseif strcmp(list1{index1}, '2022') == 1
    cd png2022
else
    cd png2023
end

% Month
d = dir;
list2 = {d(3:end).name};
[index2, dummy2] = listdlg('PromptString', {'Select a Month.'}, 'SelectionMode', 'single', 'ListString', list2);
if dummy2 == 0
    cd ../
    return;
elseif strcmp(list2{index2}, '01') == 1
    cd 01
elseif strcmp(list2{index2}, '02') == 1
    cd 02
elseif strcmp(list2{index2}, '03') == 1
    cd 03
elseif strcmp(list2{index2}, '04') == 1
    cd 04
elseif strcmp(list2{index2}, '05') == 1
    cd 05
elseif strcmp(list2{index2}, '06') == 1
    cd 06
elseif strcmp(list2{index2}, '07') == 1
    cd 07
elseif strcmp(list2{index2}, '08') == 1
    cd 08
elseif strcmp(list2{index2}, '09') == 1
    cd 09
elseif strcmp(list2{index2}, '10') == 1
    cd 10
elseif strcmp(list2{index2}, '11') == 1
    cd 11
else
    cd 12
end

% File
d = dir;
list3 = {d.name};
[index3, dummy3] = listdlg('PromptString', {'Select a File.'}, 'SelectionMode', 'single', 'ListString', list3);
if dummy3 == 0
    cd ../../
    return;
end

%% Image
pnglong = imread(list3{index3});
pnglong = pnglong(1 : 512 * 1080 * 128);
pnglong = reshape(pnglong, 512, 1080, 128);
pnglong = flip(pnglong, 2);
pnglong = flip(pnglong, 3);

% Grid
modifiy_theta = pi * 5/3;

r = linspace(800, 2333, 512);
theta = linspace(0 - modifiy_theta, 2*pi - modifiy_theta, 1080);

x = r' * cos(theta);
y = r' * sin(theta);

% Zone
zone_lenght_x = 600;
zone_length_y = 600;

zx = -zone_lenght_x/2 : 3 : zone_lenght_x/2;
zy = -zone_length_y/2 : 3 : zone_length_y/2;
[Zx_sub, Zy_sub] = meshgrid(zx, zy);

% Surf zone
surf_theta1 = deg2rad(90); % box rotation
surf_theta2 = deg2rad(165); % position rotation
surf_center = 830 + zone_length_y/2; % center distance

surf_Zx_center = Zx_sub * cos(surf_theta1) - Zy_sub * sin(surf_theta1);
surf_Zy_center = Zx_sub * sin(surf_theta1) + Zy_sub * cos(surf_theta1);

surf_Zx_real = surf_Zx_center + surf_center * cos(surf_theta2);
surf_Zy_real = surf_Zy_center + surf_center * sin(surf_theta2);

temp = surf_Zx_real + surf_Zy_real * 1i;
temp_r = abs(temp);
temp_th = angle(temp);
temp_th = mod(modifiy_theta + temp_th, 2*pi);
pngsurf = zeros(201, 201, 128);
pngsurf = single(pngsurf);

for i = 1 : 128
    sub = pnglong(:, :, i);

    idx_r = floor((temp_r - 800) / 3) + 1;
    idx_th = floor(temp_th / abs(theta(1) - theta(2))) + 1;
    idx = idx_r + 512 * idx_th;
    dummy = sub(idx);
    pngsurf(:, :, i) = dummy;
end

% Wave zone
wave_theta1 = deg2rad(90); % box rotation
wave_theta2 = deg2rad(145); % position rotation
wave_center = 1150 + zone_length_y/2; % center distance

wave_Zx_center = Zx_sub * cos(wave_theta1) - Zy_sub * sin(wave_theta1);
wave_Zy_center = Zx_sub * sin(wave_theta1) + Zy_sub * cos(wave_theta1);

wave_Zx_real = wave_Zx_center + wave_center * cos(wave_theta2);
wave_Zy_real = wave_Zy_center + wave_center * sin(wave_theta2);

temp = wave_Zx_real + wave_Zy_real * 1i;
temp_r = abs(temp);
temp_th = angle(temp);
temp_th = mod(modifiy_theta + temp_th, 2*pi);
pngwave = zeros(201, 201, 128);
pngwave = single(pngwave);

for i = 1 : 128
    sub = pnglong(:, :, i);

    idx_r = floor((temp_r - 800) / 3) + 1;
    idx_th = floor(temp_th / abs(theta(1) - theta(2))) + 1;
    idx = idx_r + 512 * idx_th;
    dummy = sub(idx);
    pngwave(:, :, i) = dummy;
end


% Figure
for i = 1 : 128
    figure(1);
    colormap default
    tiledlayout(1,3);
    nexttile
    surf(x, y, pnglong(:, :, i), 'EdgeAlpha', 0);
    %title([list3{index3}(5:8), '/', list3{index3}(9:10), '/', list3{index3}(11:12), ' ', list3{index3}(14:15), ':', list3{index3}(16:17), ' (', num2str(i), '/128)'], FontSize=15);
    xlabel('[m]'); ylabel('[m]');
    axis equal;
    axis([-2500 2500 -2500 2500]);
    view(0, 90);
    set(gcf, 'Position', [50, 980-700, 1920-100, 500]);
    hold on;
    surf_xx = [surf_Zx_real(1, 1), surf_Zx_real(1, 201); surf_Zx_real(201, 1), surf_Zx_real(201, 201)];
    surf_yy = [surf_Zy_real(1, 1), surf_Zy_real(1, 201); surf_Zy_real(201, 1), surf_Zy_real(201, 201)];
    surf(surf_xx, surf_yy, [201, 201; 201, 201], 'FaceAlpha', 0, 'EdgeColor', 'red', 'LineWidth', 2);
    wave_xx = [wave_Zx_real(1, 1), wave_Zx_real(1, 201); wave_Zx_real(201, 1), wave_Zx_real(201, 201)];
    wave_yy = [wave_Zy_real(1, 1), wave_Zy_real(1, 201); wave_Zy_real(201, 1), wave_Zy_real(201, 201)];
    surf(wave_xx, wave_yy, [201, 201; 201, 201], 'FaceAlpha', 0, 'EdgeColor', 'red', 'LineWidth', 2);
    hold off
    nexttile
    surf(surf_Zx_real, surf_Zy_real, pngsurf(:, :, i), 'EdgeAlpha', 0);
    axis equal;
    view(0, 90);
    set(gca,'xtick',[],'ytick',[]);
    title('Surf Zone', FontSize=15);
    nexttile
    surf(wave_Zx_real, wave_Zy_real, pngwave(:, :, i), 'EdgeAlpha', 0);
    axis equal;
    view(0, 90);
    set(gca,'xtick',[],'ytick',[]);
    title('Wave Zone', FontSize=15);
end

% Restore directory
cd ../../