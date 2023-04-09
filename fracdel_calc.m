clc; clear;
close all

Fs = 1e7;   % частота дискретизации 
sps = 7; % число отсчетов на символ
L = 71; % длина фильтра (количество отсчетов)
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 1; % степень сглаживания

M = 1; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

numVars = 5; % Количество точек для построения параболы

snr(:,1) = (-5:0.5:45); % Вектор ОСШ

delay = 0.4;

data = [1,0,1,0,1,1,0,0,1,1,0,1,1,1,0,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1,0,1,1,1,1,0,0,1,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,1,1,1,1,0];
data = qammod(data, modOrder, 'UnitAveragePower' , true);

h(:,1) = createH(L, Ts, T, beta);

output = upsample(data, sps);
    
signal = conv(h, output);
varSignal = var(signal);

hd = fdesign.fracdelay(delay, 'n', 3);
fir = design(hd, 'lagrange', 'FilterStructure', 'farrowfd');

ard=[];
ard2=[];
quad=[];
quad2=[];
ccf=[];
np=[];
y=[];

for i = 1:length(snr)
    for j = 1:10000
        signal2 = noiseSignal(signal, fir, varSignal, snr(i));

        [ccf, lags] = xcorr(output, signal2);

        [~, ind] = max(abs(ccf));

        [del, np, y] = calcDelayMMSE(ccf, numVars, ind, sps);

        ard(i, j) = abs(del - delay);
        ard2(i, j) = abs(calcDelayPrecise(abs(ccf), ind) - delay);
        
        quad(i, j) = (del - delay)^2;
        quad2(i, j) = (calcDelayPrecise(abs(ccf), ind) - delay)^2;
    end
end

figure(1);
plot(abs(ccf));
hold on;
plot(np, y);
legend('CCF', 'Interpolation');
hold off;

errMMSE=[];
errPrecise=[];
dispers=[];
dispers2=[];

for i = 1:length(snr)
    errMMSE(i) = sum(ard(i,:))/length(ard(i,:));
    errPrecise(i) = sum(ard2(i,:))/length(ard2(i,:));
    
    dispers(i) = (sum(quad(i,:))/length(quad(i,:))) - (sum(ard(i,:))/length(ard(i,:)))^2;
    dispers2(i) = (sum(quad2(i,:))/length(quad2(i,:))) - (sum(ard2(i,:))/length(ard2(i,:)))^2;
end

figure(2);
plot(snr, errMMSE); xlabel('SNR'); ylabel('Mean error');
hold on;
plot(snr, errPrecise);
legend('MMSE', 'precise');
grid on;
hold off;

figure(3);
plot(snr, sqrt(dispers)); xlabel('SNR'); ylabel('Mean square deviation');
hold on;
plot(snr, sqrt(dispers2));
legend('MMSE', 'precise');
grid on;
hold off;

function [del, np, y] = calcDelayMMSE(ccf, numVars, ind, sps)
    res = calcRes(ind, ccf, numVars);

    np = linspace(ind - ceil(sps/3), ind + ceil(sps/3), 3000);
    y = res(1) + np*res(2) + res(3)*np.^2;
    
    [~, argM] = max(y);

    argMax = np(argM);
    del = argMax - ind;

    if del > 0
        del = 1 - del;
    else
        del = -del;
    end
end

function [signal2] = noiseSignal(signal, fir, varSignal, snr)
    signal2 = filter(fir, signal);

    varNoise = varSignal*10^(-snr/10);
    signal2 = signal2+sqrt(varNoise/2)*(randn(size(signal2))+1i*randn(size(signal2)));
end

function [res] = calcRes(ind, ccf, numVars)
    A=[];
    b=[];
    for i = 1:numVars
        A(i, 1) = 1;
        A(i, 2) = (ind - floor(numVars/2) + i - 1);
        A(i, 3) = (ind - floor(numVars/2) + i - 1)^2;
        b(i,1) = (ccf(ind - floor(numVars/2) + i - 1));
    end
    res = inv(A.'*A)*A.'*b;
end

function [del] = calcDelayPrecise(ccf, ind)
    del = ((ccf(ind-1) - ccf(ind+1))/(2*(ccf(ind-1) + ccf(ind+1) - 2*ccf(ind))));
    if del > 0
        del = 1 - del;
    else
        del = -del;
    end
end
