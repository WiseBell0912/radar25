clear; close all; clc;

%load('test_0414.mat');
load("result_1910.mat");
load('snr_y1910_0416.mat');

%%
filter = (LandE < mean(LandE) + 2 * std(LandE));

figure;
set(gcf, 'WindowState', 'maximized');
drawnow;
tiledlayout(2, 4);

nexttile([1, 4]);
hold on;
plot(Date, surf_Peakiness);


%nexttile([1, 2]);
plot(Date, wave_Peakiness);
hold off;
legend('surf', 'wave');

% nexttile([1, 2]);
% plot(Date, rs_surf_maxW);
% title('surf max W');
% 
% nexttile([1, 2]);
% plot(Date, rs_wave_maxW);
% title('wave max W');

dcm = datacursormode(gcf);
datacursormode on;
extraParam1 = Date;
extraParam2 = img_surf;
extraParam3 = img_wave;
set(dcm, 'UpdateFcn', @(obj, event_obj) cursorCallback(obj, event_obj, extraParam1, extraParam2, extraParam3));
function txt = cursorCallback(~, event_obj, extraParam1, extraParam2, extraParam3)
persistent infoTextHandle
pos = event_obj.Position; % [x, y] 값
clickedDate = datetime(pos(1), 'ConvertFrom', 'datenum') + calyears(2019) + calmonths(9) + caldays(1); % x값을 datetime으로 변환

idx = find(extraParam1 == clickedDate);
surf_img = extraParam2(:, :, idx);
wave_img = extraParam3(:, :, idx);

% 'yyyy-mm-dd HH:MM' 포맷으로 변환하여 출력
txt = {sprintf('Date = %s', datestr(clickedDate, 'yyyy-mm-dd HH:MM')), ...
    sprintf('value = %f', pos(2))}; % 추가 인수 사용
% 두 번째 타일에 이미지 업데이트
nexttile(6);
surf(surf_img, 'EdgeAlpha', 0);
title('surf');
view(0, 90);
axis equal;
axis off;
xlim([1, 201]);
ylim([1, 201]);
nexttile(7);
surf(wave_img, 'EdgeAlpha', 0);
title('wave');
view(0, 90);
axis equal;
axis off;
xlim([1, 201]);
ylim([1, 201]);
%imshow(['/Users/limhyeonjong/Documents/Personal/GraduateProject/Image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
%imshow(['C:/Users/Hyeonjong Im/Documents/새 폴더/image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);

end
