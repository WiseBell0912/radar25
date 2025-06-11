clear; clc;

% Figure
figure(6);
tiledlayout(1, 2);

% 1
file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/AIB_20191205_1820.png';
modifiy_theta = pi * 5/3;

pnglong = imread(file_path);

pnglong = pnglong(1 : 512 * 1080 * 128);
pnglong = reshape(pnglong, 512, 1080, 128);
pnglong = flip(pnglong, 2);
pnglong = flip(pnglong, 3);

r = linspace(800, 2333, 512);
theta = linspace(0 - modifiy_theta, 2*pi - modifiy_theta, 1080);

x = r' * cos(theta);
y = r' * sin(theta);

nexttile;
for i = 1 : 1
    surf(x, y, pnglong(:, :, i), 'EdgeAlpha', 0);
    axis equal;
    grid off;
    axis([-2500 2500 -2500 2500]);
    view(0, 90);
    title('(a) Clear Day Radar Image 2019/12/05 18:20',  'FontName', 'Times New Roman', 'FontSize', 20);
    xlabel('[m]'); ylabel('[m]');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
end

% 2
file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/AIB_20200907_0820.png';
modifiy_theta = pi * 5/3;

pnglong = imread(file_path);

pnglong = pnglong(1 : 512 * 1080 * 128);
pnglong = reshape(pnglong, 512, 1080, 128);
pnglong = flip(pnglong, 2);
pnglong = flip(pnglong, 3);

r = linspace(800, 2333, 512);
theta = linspace(0 - modifiy_theta, 2*pi - modifiy_theta, 1080);

x = r' * cos(theta);
y = r' * sin(theta);

nexttile;
for i = 1 : 1
    surf(x, y, pnglong(:, :, i), 'EdgeAlpha', 0);
    axis equal;
    grid off;
    axis([-2500 2500 -2500 2500]);
    view(0, 90);
    title('(b) Rainy Day Radar Image 2020/09/07 08:20',  'FontName', 'Times New Roman', 'FontSize', 20);
    xlabel('[m]'); ylabel('[m]');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
end