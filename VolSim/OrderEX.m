clear;
%figHandle = openfig('1114_6_1.fig', 'invisible');
%figHandle = openfig('1118_5_2.fig', 'invisible');
%figHandle = openfig('1118_6_3.fig', 'invisible');
%figHandle = openfig('matlab_figure_20241001_1015_02.fig', 'invisible');
%figHandle = openfig('matlab_figure_20241001_1017_02.fig', 'invisible');
figHandle = openfig('./exports/6-1sim.fig', 'invisible');



axesHandles = findobj(figHandle, 'Type', 'axes');
lineHandles = findobj(axesHandles, 'Type', 'line');
commonX = get(lineHandles(1), 'XData')';  % 列ベクトルに変換
numLines = length(lineHandles);
numPoints = length(commonX);
yMatrix = zeros(numPoints, numLines);  % 各列に1つの時系列データが入る

for k = 1:numLines
    % 各 line のデータ取得
    xdata = get(lineHandles(k), 'XData')';
    ydata = get(lineHandles(k), 'YData')';

    
    yMatrix(:, k) = ydata;
end

order = zeros(numPoints,6);

thetas = 0:2*pi/numLines:2*pi-2*pi/numLines;
P = zeros(2,numLines);
for i = 1:numLines
    P(1,i) = sin(thetas(i));
    P(2,i) = -cos(thetas(i));
end

for i = 1:numPoints
    for j = 1:numLines
        order(i,1) = order(i,1) + P(1,j)*sin(yMatrix(i,j))/numLines;
        order(i,2) = order(i,2) + P(1,j)*cos(yMatrix(i,j))/numLines;
        order(i,4) = order(i,4) + P(2,j)*sin(yMatrix(i,j))/numLines;
        order(i,5) = order(i,5) + P(2,j)*cos(yMatrix(i,j))/numLines;
    end
    order(i,3) = sqrt(order(i,1)^2+order(i,2)^2);
    order(i,6) = sqrt(order(i,4)^2+order(i,5)^2);
end

figure;
plot(commonX,order(:,3),commonX,order(:,6))
ylim([0 1]);
xlim([0 commonX(end)])
ylabel("R")
xlabel("Time t[s]")
legend("$$R_1$$","$$R_2$$")
grid on
tuneFigure;
% tuneFigure の後にラインの太さを再設定する
lineHandles = findobj(gca, 'Type', 'line');
set(lineHandles, 'LineWidth', 5);
set(gca, 'FontSize', 24);
% 9. 作業が終了したら figure を閉じる
close(figHandle);
