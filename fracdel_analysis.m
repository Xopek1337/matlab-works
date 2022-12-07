clc; clear;

L = 40;     % длина фильтра (количество отсчетов) - любое значение, кратное
            % 20 и не меньше 40
T = 1e-6;   % длительность символа
Fs = 1e7;   % частота дискретизации 
sps = Fs*T; % число отсчетов на символ
Ts = 1/Fs;  % период дискретизации
beta = 0.9; % степень сглаживания

M = 4; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

snr(:,1) = [-10:0.5:25]; % Вектор ОСШ
eyeSNR = 30; % ОСШ для построения глазковых диаграмм

data = randi([0 1], 5000, M); % Случайная последовательность бит

% BER без дробной задержки
ber = [];

% BER с дробной задержкой
berFir = [];

% Создание фильтра дробной задержки
hd = fdesign.fracdelay(0.5, 'n', 3);
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
plot(snr, log(ber));
hold on; 
grid on; axis('tight'); 
xlabel('SNR'); ylabel('BER'); title('Without Fractional Delay');
hold off;

figure(2); 
plot(snr, log(berFir));
hold on;
grid on; axis('tight'); 
xlabel('SNR'); ylabel('BER'); title('With Fractional Delay');
hold off;

% Построение глазковых диаграмм с BPSK-модулированным сигналом,
% сглаженным фильтром Найквиста
dataSym = bi2de(data);
data2 = qammod(dataSym, modOrder, 'UnitAveragePower' , true);

output2 = upsample(data2, sps);
output2 = awgn(output2, eyeSNR, 'measured');

data3 = filter(fir, output2);

y2 = conv(h, output2);
y3 = conv(h, data3);

eyediagram(y2, 2*sps); title('Without Fractional Delay'); % Без дробной задержки
eyediagram(y3, 2*sps); title('With Fractional Delay'); % С дробной задержкой


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
    demodData = downsample(demodData, sps);
    
    dataOut = de2bi(demodData, M);
    
    % При свертке с импульсной характеристикой добавляются лишние биты в
    % начало и конец, убираем их
    rightDataOut=[];
    
    numExtraBits = L/20;
    
    for j = 1:M
        for i = (numExtraBits + 1):(length(dataOut) - numExtraBits)
            rightDataOut(i - numExtraBits, j) = dataOut(i, j);
        end
    end
    
    [nErrors, ber] = biterr(data, rightDataOut);
end
