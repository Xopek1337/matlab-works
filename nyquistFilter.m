Fs = 1e7;   % частота дискретизации 
sps = 10; % число отсчетов на символ
L = sps*3;      % длина фильтра (количество отсчетов)
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 0.9; % степень сглаживания

bits=[0,1,0,1,1,0,1,0,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,1,1,1,0,0,1,0,1,1,1,0,0,1,0,1,0,1,1,1]; % последовательность бит
data=2*bits-1; % создание BPSK символов  
nsym=length(data); % количество символов

h = createH(L,Ts,T,beta);

output=upsample(data,sps); % увеличение частоты дискретизации
y=conv(h,output); % свертка символов с ИХ фильтра

figure(1); 
subplot(3,2,[1 2]); stem(h); grid on; axis('tight'); 
title('ИХ фильтра'); ylabel('амплитуда'); xlabel('время');
subplot(3,2,[3 4]); stem(data); grid on; axis('tight');  
title('Символы на входе фильтра'); ylabel('амплитуда'); xlabel('символы'); 
subplot(3,2,[5 6]); plot(y); grid on; axis('tight');  
title('Символы на выходе фильтра');ylabel('амплитуда'); xlabel('выборки');
