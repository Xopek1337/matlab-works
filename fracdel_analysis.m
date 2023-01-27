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
delays(:,1) = [0:0.025:0.5]; % Вектор задержек для построения сигнальных созвездий
constSNR = 1000; % ОСШ для построения сигнальных созвездий

snr(:,1) = [-2:0.5:30]; % Вектор ОСШ
data = randi([0 1], 100000, M); % Случайная последовательность бит

% BER без дробной задержки
ber = [];

% BER с дробной задержкой
berFir = [];

% Создание фильтра дробной задержки
hd = fdesign.fracdelay(0.3, 'n', 3);
fir = design(hd, 'lagrange');

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

% Рассчет BER без дробной задержки
for n = 1:length(snr)
    ber(n,1) = calcBER(snr(n), data, modOrder, 0, h, sps, M, L);
end

% Рассчет BER с дробной задержкой
for n = 1:length(snr)
    berFir(n,1) = calcBER(snr(n), data, modOrder, fir, h, sps, M, L);
end

figure(1); 
semilogy(snr, ber);
hold on; 
grid on; axis('tight'); 
xlabel('SNR'); ylabel('BER'); title('Without Fractional Delay');
hold off;

figure(2); 
semilogy(snr, berFir);
hold on;
grid on; axis('tight'); 
xlabel('SNR'); ylabel('BER'); title('With Fractional Delay');
hold off;

% Построение сигнальных созвездий
for n = 1:length(delays)
    plotConstellDiag(dataConstell, delays(n), modOrder, sps, constSNR, h);
end

function [ber] = calcBER(snr, data, modOrder, fir, h, sps, M, L)
    dataSym = bi2de(data);
    modData = qammod(dataSym, modOrder, 'UnitAveragePower' , true);
  
    output = upsample(modData, sps);
    
    signal = conv(h, output);
    signal = awgn(signal, snr, 'measured');
    
    if(fir ~= 0)
        signal = filter(fir, signal);
    end
    
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
    dataSym = bi2de(data);
    modData = qammod(dataSym, modOrder, 'UnitAveragePower' , true);

    output = upsample(modData, sps);
    
    signal = conv(h, output);
    signal = awgn(signal, constSNR, 'measured');
    
    hd = fdesign.fracdelay(delay, 'n', 3);
    fir = design(hd, 'lagrange', 'FilterStructure', 'farrowfd');

    delSignal = filter(fir, signal);

    fig = scatterplot(delSignal, sps);
    
    ax = fig.CurrentAxes;
    title(ax, 'Delay: ', delay);
end
