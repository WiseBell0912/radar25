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

% Image
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

% Store frames for GIF
frames = cell(1, 128);

% Figure
for i = 1 : 128
    figure(1);
    surf(x, y, pnglong(:, :, i), 'EdgeAlpha', 0);
    title([list3{index3}(5:8), '/', list3{index3}(9:10), '/', list3{index3}(11:12), ' ', list3{index3}(14:15), ':', list3{index3}(16:17), ' (', num2str(i), '/128)'], 'FontSize', 15);
    xlabel('[m]'); ylabel('[m]');
    axis equal;
    axis([-2500 2500 -2500 2500]);
    view(0, 90);
    set(gcf, 'Position', [0, 0, 900, 900]);
    
    % Capture the plot as an image
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    frames{i} = {imind, cm};
    
    % Close the figure
    close(figure(1));
end

% Restore directory
cd ../../

% GIF File Setup
gif_filename = ['GIF', list3{index3}(4:17), '.gif'];

% Write frames to the GIF File
for i = 1:128
    [imind, cm] = frames{i}{:};
    if i == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end
