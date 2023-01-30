clc; clear;

L = 40;     % длина фильтра (количество отсчетов) - не меньше 40
Fs = 1e7;   % частота дискретизации 
sps = 10; % число отсчетов на символ
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 0.9; % степень сглаживания

M = 4; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

dataConstell = randi([0 1], 500, M); % Случайная последовательность бит для построения сигнальных созвездий
constSNR = 1000; % ОСШ для построения сигнальных созвездий

delays(:,1) = [0:0.15:0.75]; % Вектор задержек

snr(:,1) = [-2:0.25:25]; % Вектор ОСШ
data = randi([0 1], 50000, M); % Случайная последовательность бит

% BER
ber = [];

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

% Рассчет BER с дробной задержкой
for i = 1:length(delays)
    for n = 1:length(snr)
        ber(n, i) = calcBER(snr(n), data, modOrder, h, sps, M, L, delays(i));
    end
end

legendStrings = "Delay = " + string(delays);

% Построение графиков BER для разных величин дробной задержки
for i = 1:length(delays)
    semilogy(snr, ber(:, i));
    hold on; 
    grid on; axis('tight'); 
    xlabel('SNR'); ylabel('BER'); legend(legendStrings);
end
hold off;

% Построение сигнальных созвездий
for n = 1:length(delays)
    plotConstellDiag(dataConstell, delays(n), modOrder, sps, constSNR, h);
end

function [ber] = calcBER(snr, data, modOrder, h, sps, M, L, delay)
    signal = createSignal(snr, data, modOrder, h, delay, sps);
    
    demodData = qamdemod(signal, modOrder, 'UnitAveragePower' , true);
    
    rightDataOut=[];
    
    numExtraSamples = ceil(L/2);
    for i = (numExtraSamples + 1):(length(demodData) - numExtraSamples)
        rightDataOut(i -  numExtraSamples) = demodData(i);
    end
    
    demodData = downsample(rightDataOut, sps);

    dataOut = de2bi(demodData, M);
    
    [nErrors, ber] = biterr(data, dataOut);
end

function [ber] = plotConstellDiag(data, delay, modOrder, sps, constSNR, h)
    delSignal = createSignal(constSNR, data, modOrder, h, delay, sps);

    fig = scatterplot(delSignal, sps);
    
    ax = fig.CurrentAxes;
    title(ax, 'Delay: ', delay);
end
