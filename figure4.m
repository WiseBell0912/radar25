clear; close all; clc;

load("ADCP.mat");
load("BOUY.mat");

% 공통 날짜 추출
common_dates = intersect(a_Date02, b_Date);

[~, idx_a] = ismember(common_dates, a_Date02);
[~, idx_b] = ismember(common_dates, b_Date);

m_Date = common_dates;

m_a_Hs = a_Hs02(idx_a);
m_b_Hs = b_SignificantWaveHeight(idx_b);

% NaN 제거
valid_idx = ~isnan(m_a_Hs) & ~isnan(m_b_Hs);
m_a_Hs = m_a_Hs(valid_idx);
m_b_Hs = m_b_Hs(valid_idx);

% 결정계수 계산
R = corrcoef(m_a_Hs, m_b_Hs);
Rsq = R(1,2)^2;

% 시각화
figure;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

hold on;
scatter(m_a_Hs, m_b_Hs, 10, 'filled');
plot([0, 7], [0, 7], 'r--', 'LineWidth', 1)
hold off;

xlabel('ADCP', 'FontName', 'Times New Roman');
ylabel('Bouy', 'FontName', 'Times New Roman');
title(sprintf('ADCP vs BOUY (R^2 = %.3f)', Rsq), ...
      'FontName', 'Times New Roman', 'FontSize', 20);
grid on;
axis equal;
xlim([0, 7]); ylim([0, 7]);