L=40;      % длина фильтра (количество отсчетов)
R=1e6;     % скорость передачи информации (симв/с)
sps=10;     % число отсчетов на символ
Fs=sps*R;  % частота дискретизации 
T=1/R;     % длительность символа 
Ts=1/Fs;   % период дискретизации
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

function [h] = createH(L, Ts, T, beta) % фильтр Найквиста типа 'Приподнятый косинус'
    h=zeros(1,L); 
    if mod(L,2)==0     
        M=L/2;
    else
        M=(L-1)/2;
    end
    
    for n=-M:M     
        num=sin(pi*n*Ts/T)*cos(beta*pi*n*Ts/T);     
        den=(pi*n*Ts/T)*(1-(2*beta*n*Ts/T)^2);     
        h(n+M+1)=num/den;     
        if (1-(2*beta*n*Ts/T)^2)==0         
            h(n+M+1)=pi/4*sin(pi*n*Ts/T)/(pi*n*Ts/T);     
        end
        if n==0         
            h(n+M+1)=cos(beta*pi*n*Ts/T)/(1-(2*beta*n*Ts/T)^2); % sinc == 1    
        end
    end
end