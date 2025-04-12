clear; close all; clc;

%% Radar data
load("snr_y1910_0411_2.mat");
r_WaveE = WaveE;
r_wave_Uy = wave_Uy;
r_wave_Ux = wave_Ux;
r_wave_Tp = wave_Tp;
r_wave_SNR = wave_SNR;
r_wave_Pdir = wave_Pdir;
r_wave_SpectrumMax = wave_SpectrumMax;
r_SurfE = SurfE;
r_surf_Uy = surf_Uy;
r_surf_Ux = surf_Ux;
r_surf_Tp = surf_Tp;
r_surf_SNR = surf_SNR;
r_surf_Pdir = surf_Pdir;
r_surf_SpectrumMax = surf_SpectrumMax;
r_LandE = LandE;
r_Date = Date;

load("snr_y1911_0411_2.mat");
r_WaveE = [r_WaveE; WaveE];
r_wave_Uy = [r_wave_Uy; wave_Uy];
r_wave_Ux = [r_wave_Ux; wave_Ux];
r_wave_Tp = [r_wave_Tp; wave_Tp];
r_wave_SNR = [r_wave_SNR; wave_SNR];
r_wave_Pdir = [r_wave_Pdir; wave_Pdir];
r_wave_SpectrumMax = [r_wave_SpectrumMax; wave_SpectrumMax];
r_SurfE = [r_SurfE; SurfE];
r_surf_Uy = [r_surf_Uy; surf_Uy];
r_surf_Ux = [r_surf_Ux; surf_Ux];
r_surf_Tp = [r_surf_Tp; surf_Tp];
r_surf_SNR = [r_surf_SNR; surf_SNR];
r_surf_Pdir = [r_surf_Pdir; surf_Pdir];
r_surf_SpectrumMax = [r_surf_SpectrumMax; surf_SpectrumMax];
r_LandE = [r_LandE; LandE];
r_Date = [r_Date; Date];

load("snr_y1912_0411_2.mat");
r_WaveE = [r_WaveE; WaveE];
r_wave_Uy = [r_wave_Uy; wave_Uy];
r_wave_Ux = [r_wave_Ux; wave_Ux];
r_wave_Tp = [r_wave_Tp; wave_Tp];
r_wave_SNR = [r_wave_SNR; wave_SNR];
r_wave_Pdir = [r_wave_Pdir; wave_Pdir];
r_wave_SpectrumMax = [r_wave_SpectrumMax; wave_SpectrumMax];
r_SurfE = [r_SurfE; SurfE];
r_surf_Uy = [r_surf_Uy; surf_Uy];
r_surf_Ux = [r_surf_Ux; surf_Ux];
r_surf_Tp = [r_surf_Tp; surf_Tp];
r_surf_SNR = [r_surf_SNR; surf_SNR];
r_surf_Pdir = [r_surf_Pdir; surf_Pdir];
r_surf_SpectrumMax = [r_surf_SpectrumMax; surf_SpectrumMax];
r_LandE = [r_LandE; LandE];
r_Date = [r_Date; Date];

r_wave_SNR = sqrt(r_wave_SNR);
r_surf_SNR = sqrt(r_surf_SNR);

r_wave_U = sqrt(r_wave_Ux.^2 + r_wave_Uy.^2);
r_surf_U = sqrt(r_surf_Ux.^2 + r_surf_Uy.^2);

clear Date LandE surf_Pdir surf_SNR surf_SpectrumMax surf_Tp surf_Ux surf_Uy SurfE wave_Pdir wave_SNR wave_SpectrumMax wave_Tp wave_Ux wave_Uy WaveE

%% Figure
figure(1);
x_vals = linspace(min(r_surf_SpectrumMax), max(r_surf_SpectrumMax), 10000);
mu = mean(r_surf_SpectrumMax);
sigma = std(r_surf_SpectrumMax);
pdf_vals = normpdf(x_vals, mu, sigma);
plot(x_vals, pdf_vals, 'LineWidth', 2);
title('정규 분포 (PDF)');
xlabel('x'); ylabel('확률 밀도');
hold on;

% 평균선
xline(mu, '--k', 'Mean (μ)', 'LabelHorizontalAlignment','left');

% 표준편차 구간선
xline(mu - sigma, ':r', '-1σ');
xline(mu + sigma, ':r', '+1σ');
xline(mu - 2*sigma, ':g', '-2σ');
xline(mu + 2*sigma, ':g', '+2σ');
xline(mu - 3*sigma, ':b', '-3σ');
xline(mu + 3*sigma, ':b', '+3σ');

legend('PDF', 'Mean', '±1σ', '±2σ', '±3σ');
hold off;

%% Figure
figure(2);
t = tiledlayout(2,3);
a = nexttile([1 3]);
hold on;
%plot(r_Date, r_wave_SpectrumMax, 'Color', [1, 0, 0, 0.8], 'LineStyle', '-');
plot(r_Date(r_surf_SpectrumMax < mu - 0.3 * sigma), r_surf_SpectrumMax(r_surf_SpectrumMax < mu - 0.3 * sigma), 'Color', [0, 0, 1, 0.8], 'LineStyle', '-');
%yline(mu, ':r', 'Mean (μ)');
%yline(mu - 0.3 * sigma, ':r', 'Mean (μ)');
ylabel("U [m/s]");
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([datetime(2019, 10, 1), datetime(2019, 10, 31)]);
title("Current", "FontSize", 15);
xlabel("Date [mm dd]");
%legend('Wave', 'Surf');



% 데이터 커서 모드 활성화
dcm = datacursormode(gcf);
datacursormode on;

% 추가 인수 (예: 추가 정보를 담은 변수)
extraParam1 = r_Date;
extraParam2 = r_wave_U;
extraParam3 = r_surf_U;
extraParam4 = r_wave_SpectrumMax;
extraParam5 = r_surf_SpectrumMax;

% 익명 함수로 추가 인수 전달
set(dcm, 'UpdateFcn', @(obj, event_obj) cursorCallback(obj, event_obj, extraParam1, extraParam2, extraParam3, extraParam4, extraParam5));

b = nexttile([1 3]);

% 클릭 시 좌표를 출력하는 함수
function txt = cursorCallback(~, event_obj, extraParam1, extraParam2, extraParam3, extraParam4, extraParam5)
persistent infoTextHandle
pos = event_obj.Position; % [x, y] 값
clickedDate = datetime(pos(1), 'ConvertFrom', 'datenum') + calyears(2019) + calmonths(9) + caldays(1); % x값을 datetime으로 변환

idx = extraParam1 == clickedDate;
fr_wave_U = extraParam2(idx);
fr_surf_U = extraParam3(idx);
fr_wave_SpectrumMax = extraParam4(idx);
fr_surf_SpectrumMax = extraParam4(idx);


% 'yyyy-mm-dd HH:MM' 포맷으로 변환하여 출력
txt = {sprintf('Date = %s', datestr(clickedDate, 'yyyy-mm-dd HH:MM')), ...
    sprintf('U = %.2f [m/s]', pos(2))}; % 추가 인수 사용

% 두 번째 타일에 이미지 업데이트
nexttile(4, [1, 2]);
%imshow(['/Users/limhyeonjong/Documents/Personal/GraduateProject/Image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
imshow(['C:/Users/Hyeonjong Im/Documents/새 폴더/image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
nexttile(6);
axis off; % 축 숨기기
information = {...
    sprintf('Wave U  =  %6.2f  [m/s]', fr_wave_U), ...
    sprintf('Surf U  =  %6.2f  [m/s]\n', fr_surf_U), ...
    sprintf('Wave Max  =  %6.8f  []', fr_wave_SpectrumMax), ...
    sprintf('Surf Max  =  %6.8f  []', fr_surf_SpectrumMax), ...
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

