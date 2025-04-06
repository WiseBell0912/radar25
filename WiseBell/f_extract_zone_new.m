function [pngsurf, pngwave, pngenergy] = f_extract_zone(pnglong, modifiy_theta, surf_theta1, surf_theta2, surf_center, wave_theta1, wave_theta2, wave_center, energy_theta1, energy_theta2, energy_center)

% Grid
r = linspace(800, 2333, 512);
theta = linspace(0 - modifiy_theta, 2*pi - modifiy_theta, 1080);

% Zone
zone_lenght_x = 600;
zone_length_y = 600;

zx = -zone_lenght_x/2 : 3 : zone_lenght_x/2;
zy = -zone_length_y/2 : 3 : zone_length_y/2;
[Zx_sub, Zy_sub] = meshgrid(zx, zy);

% Surf zone
surf_Zx_center = Zx_sub * cos(surf_theta1) - Zy_sub * sin(surf_theta1);
surf_Zy_center = Zx_sub * sin(surf_theta1) + Zy_sub * cos(surf_theta1);

surf_Zx_real = surf_Zx_center + surf_center * cos(surf_theta2);
surf_Zy_real = surf_Zy_center + surf_center * sin(surf_theta2);

temp = surf_Zx_real + surf_Zy_real * 1i;
temp_r = abs(temp);
temp_th = angle(temp);
temp_th = mod(modifiy_theta + temp_th, 2*pi);
pngsurf = zeros(201, 201, 128);
pngsurf = double(pngsurf);

for i = 1 : 128
    sub = pnglong(:, :, i);

    idx_r = floor((temp_r - 800) / 3) + 1;
    idx_th = floor(temp_th / abs(theta(1) - theta(2))) + 1;
    idx = idx_r + 512 * idx_th;
    dummy = sub(idx);
    pngsurf(:, :, i) = dummy;
end

% Wave zone
wave_Zx_center = Zx_sub * cos(wave_theta1) - Zy_sub * sin(wave_theta1);
wave_Zy_center = Zx_sub * sin(wave_theta1) + Zy_sub * cos(wave_theta1);

wave_Zx_real = wave_Zx_center + wave_center * cos(wave_theta2);
wave_Zy_real = wave_Zy_center + wave_center * sin(wave_theta2);

temp = wave_Zx_real + wave_Zy_real * 1i;
temp_r = abs(temp);
temp_th = angle(temp);
temp_th = mod(modifiy_theta + temp_th, 2*pi);
pngwave = zeros(201, 201, 128);
pngwave = double(pngwave);

for i = 1 : 128
    sub = pnglong(:, :, i);

    idx_r = floor((temp_r - 800) / 3) + 1;
    idx_th = floor(temp_th / abs(theta(1) - theta(2))) + 1;
    idx = idx_r + 512 * idx_th;
    dummy = sub(idx);
    pngwave(:, :, i) = dummy;
end

% Energy zone
% resetting for energy zone
zone_lenght_x = 200;
zone_length_y = 500;

zx = -zone_lenght_x/2 : 3 : zone_lenght_x/2;
zy = -zone_length_y/2 : 3 : zone_length_y/2;
[Zx_sub, Zy_sub] = meshgrid(zx, zy);

energy_Zx_center = Zx_sub * cos(energy_theta1) - Zy_sub * sin(energy_theta1);
energy_Zy_center = Zx_sub * sin(energy_theta1) + Zy_sub * cos(energy_theta1);

energy_Zx_real = energy_Zx_center + energy_center * cos(energy_theta2);
energy_Zy_real = energy_Zy_center + energy_center * sin(energy_theta2);

temp = energy_Zx_real + energy_Zy_real * 1i;
temp_r = abs(temp);
temp_th = angle(temp);
temp_th = mod(modifiy_theta + temp_th, 2*pi);
pngenergy = zeros(length(zy), length(zx), 128);
pngenergy = single(pngenergy);

for i = 1 : 128
    sub = pnglong(:, :, i);

    idx_r = floor((temp_r - 800) / 3) + 1;
    idx_th = floor(temp_th / abs(theta(1) - theta(2))) + 1;
    idx = idx_r + 512 * idx_th;
    dummy = sub(idx);
    pngenergy(:, :, i) = dummy;
end

end