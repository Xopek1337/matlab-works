clc; clear;
close all

Fs = 1e7;   % частота дискретизации 
sps = 10; % число отсчетов на символ
L = 71; % длина фильтра (количество отсчетов)
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 1; % степень сглаживания

M = 4; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

dataConstell = randi([0 1], 500, M); % Случайная последовательность бит для построения сигнальных созвездий
constSNR = 1000; % ОСШ для построения сигнальных созвездий

delays(:,1) = (0.0:0.1:0.9); % Вектор задержек
snr(:,1) = (0:1:15); % Вектор ОСШ

nBits = 10000;
nRealiz = 20;

% BER
nErr = zeros(length(snr),length(delays),nRealiz);

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

% Рассчет BER с дробной задержкой
for i = 1:length(delays)
    for n = 1:length(snr)
        for iRealiz = 1:nRealiz
            data = randi([0 1], nBits, M); % Случайная последовательность бит
            nErr(n, i, iRealiz) = calcBER(snr(n), data, modOrder, h, sps, M, L, delays(i));
        end
    end
end
ber = sum(nErr,3)./(nRealiz*nBits);

legendStrings = "Delay = " + string(delays);

% Построение графиков BER для разных величин дробной задержки
figure(1);
for i = 1:length(delays)
    semilogy(snr, ber(:, i));
    hold on; 
    grid on; axis('tight'); 
    xlabel('SNR'); ylabel('BER'); legend(legendStrings);
end
hold off;

% Построение сигнальных созвездий
for n = 1:length(delays)
    figure(n+1);
    plotConstellDiag(dataConstell, delays(n), modOrder, sps, constSNR, h, L);
end

function [nErrors] = calcBER(snr, data, modOrder, h, sps, M, L, delay)
    data(1:7) = [1 0 1 0 1 1 1];
    signal = createSignal(snr, data, modOrder, h, delay, sps);

    rightDataOut = signal( (L+1)/2+1:end );

    signal = downsample(rightDataOut, sps);

    demodData = qamdemod(signal(1:length(data)), modOrder, 'UnitAveragePower' , true);
    
    dataOut = de2bi(demodData, M);
    
    [nErrors, ~] = biterr(data, dataOut);

    if snr==250
        figure;
        plot(signal,'.'); grid on
    end
    1;
end

function [ber] = plotConstellDiag(data, delay, modOrder, sps, constSNR, h, L)
    delSignal = createSignal(constSNR, data, modOrder, h, delay, sps);
    
    rightDataOut = delSignal( (L+1)/2+1:end );
    
    signal = downsample(rightDataOut, sps);
    
    plot(signal, '.');
    
    title('Delay: ', delay);
end
