clear; close all; clc;

d = dir;

year = {'2018', '2019', '2020', '2021', '2022', '2023'};
month = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'};

for i = 15 : length(d)
    disp([i, length(d)]);
    for j = 1 : 12
        if strcmp(d(i).name(9:10), month(j)) == 1
            movefile(d(i).name, ['D:\png', year{6}, '\', month{j}]);
        end
    end
end