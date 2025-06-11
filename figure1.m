clear; close all; clc;

g = 9.81;
h = 20;

kx = -2 : 0.01 : 2;
ky = -2 : 0.01 : 2;
[Kx, Ky] = meshgrid(kx, ky);
K = sqrt(Kx.^2 + Ky.^2);

W = sqrt(g .* K .* tanh(K .* h));
W2 = W + Kx .* 1;

figure(1);
tiledlayout(1, 2);

nexttile;
surf(Kx, Ky, W, 'EdgeAlpha', 0.1);
axis equal;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);  % 현재 axes 전체 텍스트에 적용
xlabel('Kx'); ylabel('Ky'); zlabel('W');
title('(a) Dispersion shell without current');
xlim([-2.5, 2.5]); ylim([-2.5, 2.5]);
zlim([min(W2(:)), max(W2(:))]);

nexttile;
surf(Kx, Ky, W2, 'EdgeAlpha', 0.1);
axis equal;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);  % 현재 axes 전체 텍스트에 적용
xlabel('Kx'); ylabel('Ky'); zlabel('W');
title('(b) Dispersion shell with X-axis current');
xlim([-2.5, 2.5]); ylim([-2.5, 2.5]);
zlim([min(W2(:)), max(W2(:))]);
