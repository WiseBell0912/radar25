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
load("BOUY.mat");
b_Hs = b_SignificantWaveHeight;
b_Wind = b_WindSpeed;
b_Pdir = b_WaveDirection;
b_Tp = b_MaximumWavePeriod;

b_Date = b_Date(~isnan(b_Hs));
b_Hs = b_Hs(~isnan(b_Hs));
b_Wind = b_Wind(~isnan(b_Hs));
b_Pdir = b_Pdir(~isnan(b_Hs));
b_Tp = b_Tp(~isnan(b_Hs));

b_Hs(b_Hs == 0) = NaN;
b_Hs = fillmissing(b_Hs, 'linear');

% ASOS
load("ASOS.mat");

% Radar wave
load("snr_y1910_0406.mat");
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

load("snr_y1911_0406.mat");
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

load("snr_y1912_0406.mat");
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
mask1 = ismember(r_Date, b_Date);
mask2 = ismember(r_Date, a_Date);
mask3 = ismember(r_Date, asos_Date);
mask = mask1 & mask2 & mask3;

rr_Date = r_Date(mask);
rr_wave_Pdir = r_wave_Pdir(mask);
rr_wave_SNR = r_wave_SNR(mask);
rr_wave_Tp = r_wave_Tp(mask);
rr_surf_Pdir = r_surf_Pdir(mask);
rr_surf_SNR = r_surf_SNR(mask);
rr_surf_Tp = r_surf_Tp(mask);

mask1 = ismember(b_Date, r_Date);
mask2 = ismember(b_Date, a_Date);
mask3 = ismember(b_Date, asos_Date);
mask = mask1 & mask2 & mask3;

bb_Date = b_Date(mask);
bb_Hs = b_Hs(mask);
bb_Wind = b_Wind(mask);
bb_Pdir = b_Pdir(mask);
bb_Tp = b_Tp(mask);

mask1 = ismember(a_Date, r_Date);
mask2 = ismember(a_Date, b_Date);
mask3 = ismember(a_Date, asos_Date);
mask = mask1 & mask2 & mask3;

aa_Date = a_Date(mask);
aa_Pdir = a_Pdir(mask);
aa_Tp = a_Tp(mask);

mask1 = ismember(asos_Date, r_Date);
mask2 = ismember(asos_Date, b_Date);
mask3 = ismember(asos_Date, a_Date);
mask = mask1 & mask2 & mask3;

aasos_Date = asos_Date(mask);
aasos_Precipitation = asos_Precipitation(mask);
%% 개발
modelfun = @(x, r2b_wave_SNR) x(1) + x(2) * r2b_wave_SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
x = lsqcurvefit(modelfun, initial_guess, rr_wave_SNR, bb_Hs, [], [], options);

rr_wave_Hs = x(1) + x(2) .* rr_wave_SNR;

modelfun = @(y, r2b_surf_SNR) y(1) + y(2) * r2b_surf_SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
y = lsqcurvefit(modelfun, initial_guess, rr_surf_SNR, bb_Hs, [], [], options);

rr_surf_Hs = y(1) + y(2) .* rr_surf_SNR;

%% 확인 Graph
figure(1);
t = tiledlayout(2,3);
a = nexttile([1 3]);
hold on;
% yyaxis left
plot(bb_Date, movmean(bb_Hs, 1), 'Color', [0, 0, 0, 0.8], 'LineStyle', '-');
plot(rr_Date, movmean(rr_wave_Hs, 1), 'Color', [1, 0, 0, 0.8], 'LineStyle', '-');
plot(rr_Date, movmean(rr_surf_Hs, 1), 'Color', [0, 0, 1, 0.8], 'LineStyle', '-');
ylabel("Hs [m]");
% yyaxis right
% plot(b2r_Date, movmean(b2r_Wind, 3), 'Color', [0.9, 0.2, 0.9, 0.7], 'LineStyle', '-.');
% ylabel("Wind Velocity [m/s]");
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([datetime(2019, 10, 1), datetime(2019, 10, 4)]);
title("Significant Wave Height", "FontSize", 15);
xlabel("Date [mm dd]");
legend('Bouy', 'Wave', 'Surf');
% legend('Bouy', 'Wave', 'Surf', 'Wind Velocity');



% 데이터 커서 모드 활성화
dcm = datacursormode(gcf);
datacursormode on;

% 추가 인수 (예: 추가 정보를 담은 변수)
extraParam1 = bb_Date;
extraParam2 = bb_Wind;
extraParam3 = rr_wave_Pdir;
extraParam4 = rr_wave_Tp;
extraParam5 = rr_surf_Pdir;
extraParam6 = rr_surf_Tp;
extraParam7 = bb_Pdir;
extraParam8 = bb_Tp;
extraParam9 = aa_Pdir;
extraParam10 = aa_Tp;
extraParam11 = aasos_Precipitation;

% 익명 함수로 추가 인수 전달
set(dcm, 'UpdateFcn', @(obj, event_obj) cursorCallback(obj, event_obj, extraParam1, extraParam2, extraParam3, extraParam4, extraParam5, extraParam6, extraParam7, extraParam8, extraParam9, extraParam10, extraParam11));

b = nexttile([1 3]);

% 클릭 시 좌표를 출력하는 함수
function txt = cursorCallback(~, event_obj, extraParam1, extraParam2, extraParam3, extraParam4, extraParam5, extraParam6, extraParam7, extraParam8, extraParam9, extraParam10, extraParam11)
persistent infoTextHandle
pos = event_obj.Position; % [x, y] 값
clickedDate = datetime(pos(1), 'ConvertFrom', 'datenum') + calyears(2019) + calmonths(9) + caldays(1); % x값을 datetime으로 변환

idx = extraParam1 == clickedDate;
fb_Wind = extraParam2(idx);
fr_wave_Pdir = extraParam3(idx);
fr_wave_Tp = extraParam4(idx);
fr_surf_Pdir = extraParam5(idx);
fr_surf_Tp = extraParam6(idx);
fb_Pdir = extraParam7(idx);
fb_Tp = extraParam8(idx);
fa_Pdir = extraParam9(idx);
fa_Tp = extraParam10(idx);
fasos_Precipitation = extraParam11(idx);

% 'yyyy-mm-dd HH:MM' 포맷으로 변환하여 출력
txt = {sprintf('Date = %s', datestr(clickedDate, 'yyyy-mm-dd HH:MM')), ...
    sprintf('Hs = %.2f [m]', pos(2))}; % 추가 인수 사용

% 두 번째 타일에 이미지 업데이트
nexttile(4, [1, 2]);
%imshow(['/Users/limhyeonjong/Documents/Personal/GraduateProject/Image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
imshow(['C:/Users/Hyeonjong Im/Documents/새 폴더/image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
nexttile(6);
axis off; % 축 숨기기
information = {...
    sprintf('Bouy Wind  =  %6.2f  [m/s]\n', fb_Wind), ...
    sprintf('ASOS Rain  =  %7.2f  [mm]\n', fasos_Precipitation), ...
    sprintf('Bouy Pdir  =  %6.2f  [deg]', fb_Pdir), ...
    sprintf('ADCP Pdir  =  %6.2f  [deg]', fa_Pdir), ...
    sprintf('Wave Pdir  =  %6.2f  [deg]', fr_wave_Pdir), ...
    sprintf('Surf Pdir  =  %6.2f  [deg]\n', fr_surf_Pdir), ...
    sprintf('Bouy Tp  =  %6.2f  [s]', fb_Tp), ...
    sprintf('ADCP Tp  =  %6.2f  [s]', fa_Tp), ...
    sprintf('Wave Tp  =  %6.2f  [s]', fr_wave_Tp), ...
    sprintf('Surf Tp  =  %6.2f  [s]', fr_surf_Tp), ...
    };
if isempty(infoTextHandle) || ~isvalid(infoTextHandle)
    % 텍스트가 아직 없으면 새로 생성
    infoTextHandle = text(0.5, 0.5, information, ...
        'Units', 'normalized', ...
        'FontName', 'Courier New', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 16 ...
        );
else
    % 기존 텍스트 객체의 내용을 업데이트
    set(infoTextHandle, 'String', information);
end
end
