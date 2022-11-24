clc; clear;

L = 40;      % длина фильтра (количество отсчетов)
sps = 10;     % число отсчетов на символ
T = 1e-6;     % длительность символа
Fs = sps/T;  % частота дискретизации 
Ts = 1/Fs;   % период дискретизации
beta = 0.9; % степень сглаживания

bits=[0,1,0,1,1,0,1,0,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,1,1,1,0,0,1,0,1,1,1,0,0,1,0,1,0,1,1,1]; % последовательность бит
bits2 = randi([0 1], 1000, 1);

dataQ = randi([0 3], 1000, 1);
data16 = randi([0 15], 1000, 1);
data64 = randi([0 63], 1000, 1);

snr = [-5:0.5:20]; % Вектор ОСШ

% BER для разных видов модуляции без дробной задержки
berQ = [];
ber16 = [];
ber64 =[ ];

% BER для разных видов модуляции с дробной задержкой
berQFir = [];
ber16Fir = [];
ber64Fir = [];

data = 2*bits-1; % создание BPSK символов  
nsym = length(data); % количество символов

h = createH(L, Ts, T, beta);

output = upsample(data, sps); % увеличение частоты дискретизации
y = conv(h, output); % свертка символов с ИХ фильтра

% Рассчет BER без дробной задержки
for n = 1:length(snr)
    berQ = horzcat(berQ, calcBER(snr(n),dataQ,4, 0));
    ber16 = horzcat(ber16, calcBER(snr(n),data16,16, 0));
    ber64 = horzcat(ber64, calcBER(snr(n),data64,64, 0));
end

% Создание фильтра дробной задержки
hd = fdesign.fracdelay(0.2, 'n', 3);
fir = design(hd, 'lagrange');

% Рассчет BER с дробной задержкой
for n = 1:length(snr)
    berQFir = horzcat(berQFir, calcBER(snr(n), dataQ, 4, fir));
    ber16Fir = horzcat(ber16Fir, calcBER(snr(n), data16, 16, fir));
    ber64Fir = horzcat(ber64Fir, calcBER(snr(n), data64, 64, fir));
end

figure(1); 
subplot(3, 2, [1 2]); stem(h); grid on; axis('tight'); 
title('ИХ фильтра'); ylabel('амплитуда'); xlabel('время');
subplot(3, 2, [3 4]); stem(data); grid on; axis('tight');  
title('Символы на входе фильтра'); ylabel('амплитуда'); xlabel('символы'); 
subplot(3, 2, [5 6]); plot(y); grid on; axis('tight');  
title('Символы на выходе фильтра'); ylabel('амплитуда'); xlabel('выборки');

figure(2); 
plot(snr, berQ);
hold on;
plot(snr, ber16);
plot(snr, ber64); 
grid on; axis('tight'); 
legend('QPSK', '16-QAM', '64-QAM'); xlabel('SNR'); ylabel('BER');
hold off;

figure(3); 
plot(snr, berQFir);
hold on;
plot(snr, ber16Fir);
plot(snr, ber64Fir); 
grid on; axis('tight'); 
legend('QPSK', '16-QAM', '64-QAM'); xlabel('SNR'); ylabel('BER'); 
hold off;

% Построение глазковых диаграмм с BPSK-модулированным сигналом,
% сглаженным фильтром Найквиста
data2 = 2*bits2 - 1;
data2 = awgn(data2, 40);

data3 = filter(fir, data2);

output2 = upsample(data2, sps);
y2 = conv(h, output2);

output3 = upsample(data3, sps);
y3 = conv(h, output3);

eyediagram(y2, 2*sps);
eyediagram(y3, 2*sps);

function [h] = createH(L, Ts, T, beta) % фильтр Найквиста типа 'Приподнятый косинус'
    h = zeros(1, L); 
    if mod(L, 2)==0     
        M = L/2;
    else
        M = (L-1)/2;
    end
    
    for n = -M:M     
        num = sin(pi*n*Ts/T)*cos(beta*pi*n*Ts/T);     
        den = (pi*n*Ts/T)*(1-(2*beta*n*Ts/T)^2);     
        h(n+M+1) = num/den;     
        if (1-(2*beta*n*Ts/T)^2) == 0         
            h(n+M+1) = pi/4*sin(pi*n*Ts/T)/(pi*n*Ts/T);     
        end
        if n == 0         
            h(n+M+1) = cos(beta*pi*n*Ts/T)/(1-(2*beta*n*Ts/T)^2); % sinc == 1    
        end
    end
end

function [ber] = calcBER(snr, data, modOrder, fir)
    cnt = 0;
    isFiltered = 0;
    modData = qammod(data, modOrder);
    modData = awgn(modData, snr);
    
    if(fir ~= 0)
        modData = filter(fir, modData);
        isFiltered = 1; %Первое значение инициализируется нулем, поэтому в отфильтрованном сигнале начинаю счет со 2 элемента
    end
    
    demodData = qamdemod(modData, modOrder);
    
    for i = 1:length(data) - isFiltered
       if(data(i) == demodData(i + isFiltered))
           cnt = cnt + 1;
       end
    end
    ber = 1 - (cnt / length(data));
end