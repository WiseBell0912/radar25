function [] = f_extract_PPI(png_long, modifiy_theta)

% Grid
r = linspace(800, 2333, 512);
theta = linspace(0 - modifiy_theta, 2*pi - modifiy_theta, 1080);

x = r' * cos(theta);
y = r' * sin(theta);

% Figure
for i = 1 : 1
    figure(1);
    surf(x, y, png_long(:, :, i), 'EdgeAlpha', 0);
    axis equal;
    axis([-2500 2500 -2500 2500]);
    view(0, 90);
    hold on;
    surf(x(:, 365:385, i), y(:, 365:385, i), png_long(:, 365:385, i));
    hold off;
end