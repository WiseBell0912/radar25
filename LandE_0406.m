clear; close all; clc;

load("ASOS.mat");
load("snr_y1910_0406.mat");
r_LandE = LandE;
r_SurfE = SurfE;
r_WaveE = WaveE;
r_Date = Date;
load("snr_y1911_0406.mat");
r_LandE = [r_LandE; LandE];
r_SurfE = [r_SurfE; SurfE];
r_WaveE = [r_WaveE; WaveE];
r_Date = [r_Date; Date];
load("snr_y1912_0406.mat");
r_LandE = [r_LandE; LandE];
r_SurfE = [r_SurfE; SurfE];
r_WaveE = [r_WaveE; WaveE];
r_Date = [r_Date; Date];

mask = ismember(r_Date, asos_Date);
rr_Date = r_Date(mask);
rr_LandE = r_LandE(mask);
rr_SurfE = r_SurfE(mask);
rr_WaveE = r_WaveE(mask);

mask = ismember(asos_Date, r_Date);
aasos_Date = asos_Date(mask);
aasos_Precipitation = asos_Precipitation(mask);

figure(1);
t = tiledlayout(2,1);
a = nexttile([1 1]);
hold on;

% 현재 y축의 한계를 얻기 위해 임의로 플롯 (아래에서 patch 적용 전에 ylimit이 설정되어야 함)
h1 = plot(rr_Date, rr_LandE);
h2 = plot(rr_Date, rr_SurfE);
h3 = plot(rr_Date, rr_WaveE);
%h4 = plot(aasos_Date, ones(size(rr_LandE))*(4.5*10^7));
legend('land', 'surf', 'wave');

% y축 한계를 가져옵니다.
yl = ylim;

% 0이 아닌 인덱스를 찾습니다.
nonzeroIdx = find(aasos_Precipitation ~= 0);

% 연속된 인덱스 구간(group) 찾기
groups = {};
if ~isempty(nonzeroIdx)
    group = nonzeroIdx(1);
    for k = 2:length(nonzeroIdx)
        % 만약 인덱스가 연속이면 group에 추가, 아니면 그룹을 저장하고 새 그룹 시작
        if nonzeroIdx(k) == nonzeroIdx(k-1) + 1
            group = [group; nonzeroIdx(k)];
        else
            groups{end+1} = group;
            group = nonzeroIdx(k);
        end
    end
    groups{end+1} = group;  % 마지막 그룹 저장
end

% 각 그룹별로 배경 채우기 (형광 노란색, 투명도 0.3)
% 시간 간격을 고려하여 각 구간의 시작과 끝을 구합니다.
% 여기서는 인접한 날짜간 간격이 일정하다고 가정하고,
% 양쪽에 절반 간격씩 확장하여 표시합니다.
if length(aasos_Date) > 1
    dt = median(diff(aasos_Date));  % 시간 간격 (날짜 형식이면 단위: days)
else
    dt = 0;
end

for i = 1:length(groups)
    idxGroup = groups{i};
    x_start = aasos_Date(idxGroup(1)) - dt/2;
    x_end   = aasos_Date(idxGroup(end))  + dt/2;
    
    % patch를 이용해 사각형 영역 채우기
    patch([x_start, x_end, x_end, x_start], [yl(1), yl(1), yl(2), yl(2)], ...
        [0.9 0.9 0.1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% patch가 배경에 깔리도록 플롯 순서를 조정
uistack(h1, 'top');
uistack(h2, 'top');
uistack(h3, 'top');

hold off;

% 데이터 커서 모드 활성화
dcm = datacursormode(gcf);
datacursormode on;

set(dcm, 'UpdateFcn', @(obj, event_obj) cursorCallback(obj, event_obj));

% 클릭 시 좌표를 출력하는 함수
function txt = cursorCallback(~, event_obj)
pos = event_obj.Position; % [x, y] 값
clickedDate = datetime(pos(1), 'ConvertFrom', 'datenum') + calyears(2019) + calmonths(9) + caldays(1); % x값을 datetime으로 변환

% 'yyyy-mm-dd HH:MM' 포맷으로 변환하여 출력
txt = {sprintf('.')}; % 추가 인수 사용

% 두 번째 타일에 이미지 업데이트
nexttile(2);
%imshow(['/Users/limhyeonjong/Documents/Personal/GraduateProject/Image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
imshow(['C:/Users/Hyeonjong Im/Documents/새 폴더/image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
end