clc; clear;

Fs = 1e7;   % частота дискретизации 
sps = 10; % число отсчетов на символ
L = 71; % длина фильтра (количество отсчетов)
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 1; % степень сглаживания

M = 1; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

delay = 0.7;

snr = 1000;
data = randi([-1 1], 63, M);
h(:,1) = createH(L, Ts, T, beta);

output = upsample(data, sps);
    
signal = conv(h, output);
varSignal = var(signal);

hd = fdesign.fracdelay(delay, 'n', 3);
fir = design(hd, 'lagrange', 'FilterStructure', 'farrowfd');
    
signal2 = filter(fir, signal);

varNoise = varSignal*10^(-snr/10);
signal2 = signal2+sqrt(varNoise/2)*(randn(size(signal2))+1i*randn(size(signal2)) );

[ccf, lags] = xcorr(signal, signal2);

figure(1);
plot(lags, abs(ccf));

[~, ind] = max(abs(ccf));

% Сравнение вычисленной задержки с заданной
disp('calcucalted delay: ');
disp(calcDelay(abs(ccf), ind));
disp('given delay: ');
disp(delay);

function [del] = calcDelay(ccf, ind)
    del = ((ccf(ind-1) - ccf(ind+1))/(2*(ccf(ind-1) + ccf(ind+1) - 2*ccf(ind))));
    if del > 0
        del = 1 - del;
    else
        del = -del;
    end
end