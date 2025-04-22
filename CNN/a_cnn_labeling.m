%% Surf zone과 Wave zone의 파도 발생 여부를 라벨링하는 작업을 편하게 도와주는 코드.
clear; close all; clc;

%% -- 1) 라벨을 저장할 배열 초기화 --
% PNG 파일 목록 구하기
imageDir = 'C:\Users\Hyeonjong Im\Documents\새 폴더\image201910\';       % 사용자의 이미지 폴더 경로
files = dir(fullfile(imageDir, '*.png'));
nFiles = numel(files);

if nFiles == 0
    error('PNG 파일이 폴더 내에 없습니다: %s', imageDir);
end

% 라벨을 저장할 셀/문자/문자열 배열
labelCell = cell(nFiles, 1);

%% -- 2) Figure 생성 --
fig = figure('Name','Manual Labeling','Color','w',...
    'Units','normalized','Position',[0.2,0.2,0.4,0.6]);

for i = 1:nFiles

    % -- 2-1) 이미지 읽기 --
    imgPath = fullfile(files(i).folder, files(i).name);
    I = imread(imgPath);

    % (옵션) 레이더라면 grayscale(단일 채널)일 텐데,
    %        uint8 / double 등 타입이 다양할 수 있으므로
    %        imshow(I,[]) 형태로 표시
    imshow(I, []);
    title(sprintf('Image %d / %d : %s', i, nFiles, files(i).name), 'Interpreter','none');
    drawnow;

    % -- 2-2) 라벨 받기 --
    % (A) questdlg(버튼)로 라벨 받기
    choice = questdlg('Is this wave or noWave?', ...
        'Labeling', 'wave','noWave','Cancel','wave');

    % 사용자가 'Cancel' 클릭하거나 창 닫은 경우 → 중단
    if isempty(choice) || strcmpi(choice,'Cancel')
        disp('사용자가 라벨링을 취소했습니다. 현재까지 라벨 저장 후 종료합니다.');
        break;
    end

    % 라벨 기록
    labelCell{i} = choice;

    % (B) 만약 키보드로 직접 입력받고 싶다면 (A)를 주석 처리하고 사용
    %{
        prompt = sprintf('[%d/%d] Enter label (wave / noWave): ', i, nFiles);
        userLabel = input(prompt, 's');
        if isempty(userLabel)
            disp('사용자가 라벨링을 취소했습니다. 저장 후 종료합니다.');
            break;
        end
        labelCell{i} = userLabel;
    %}

end

% i-1번째까지 라벨을 달았을 수 있으므로, 실제 라벨링된 개수
labeledCount = sum(~cellfun(@isempty,labelCell));
fprintf('총 %d개 중 %d개 라벨링 완료\n', nFiles, labeledCount);

%% -- 3) 라벨 데이터 저장 --
% 3-1) MAT 파일로 저장
labelVec = categorical(labelCell(1:499));  % 분류용이므로 categorical
save('labeledPngData.mat','labelVec','-v7.3');
fprintf('labeledPngData.mat에 labelVec를 저장했습니다.\n');

close(fig);
disp('라벨링 프로세스 종료');