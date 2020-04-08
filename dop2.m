clc
clear all
close all
nfig = 1;

k = 8;
g = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1];
gsize = length(g) - 1;  %length CRC
n = k + gsize;  

mes = zeros(2^k, k);
for i = 0 : 2^k - 1
    mes(i + 1, :) = de2bi(i, k);	%Генерирование всех сообщений
end              

codes = zeros(2^k, n);


for i = 0 : 2^k - 1
    codes(i + 1, :) = [zeros(1, gsize), mes(i + 1, :)];
    [~, r] = gfdeconv(codes(i + 1, :), g);									%Остаток от деления нашего многочлена на порождающий многочлен
    codes(i + 1, :) = xor(codes(i + 1, :), [r, zeros(1, n - length(r))]);	%Кодовое слово
end


mes_by_3 = zeros(2^k, 3, 11); 
mes_modBPSK = zeros(2^k, 3, 12);
mod_Ham = zeros (2^k, 3, 15);
[H, G] = hammgen(4); %Проверочная и порождающая матрицы

for i = 0 : 2^k - 1
    mes_by_3(i + 1, :, :) = [reshape(codes(i + 1, :), 8, 3)'...
                            [0 0 0; 0 0 0; 0 0 0]];  								%Разделяем кодовое слова на 3 части и добавляем в каждую 000 в начало
    
    for j = 1 : 3
        mod_Ham(i + 1, j, :) = mod(reshape(mes_by_3(i + 1, j, :), 1, 11) * G, 2);	%Кодируем по Хэммингу
        mes_modBPSK(i + 1, j, :) = mod_Ham(i + 1, j, 1:end-3)*(-2) + 1; 			%Удаляем нули и модулируем по BPSK
    end
end                  
    

SNRdB = -10:10;
SNR = 10.^(SNRdB./10);
Pe_bit_theor = qfunc(sqrt(2.*SNR));         		%Теоретическая вероятность ошибки на бит
Ped_theor_up = ones(1, length(SNR)).*1/2^gsize; 	%Асимптотическая верхняя граница вероятности ошибки декодирования
Ped_theor = zeros(1,length(SNRdB));					%Точное значение вероятности ошибки декодирования


A = zeros(1, length(codes(1, :)) + 1);		%Книга весов
for i = 1: length(codes(:,1))
    w = sum(codes(i, :));
    A(1, w + 1) = A(1, w + 1) + 1;			
end

d = min(sum(codes(2:end, :), 2));			%Минимальное расстояние кода

for i = 1 : length(SNR)
    for j = d : n
        Ped_theor(1, i) = Ped_theor(1, i) + A(j + 1) * Pe_bit_theor(i)^j * (1 - Pe_bit_theor(i))^(n - j);
    end
end

all_codewords = de2bi((0:(2^8-1))');
h = zeros(2^8, 12);
for j=1:2^8
    temp = mod([all_codewords(j, :) 0 0 0] * G, 2);
    h(j, :) = temp(1, 1:end-3) .*(-2) + 1;			
end


Pe_bit = zeros(1,length(SNRdB));                    %Вероятность ошибки на бит
Pe_bit_after_Hamming = zeros(1,length(SNRdB));      %Вероятность ошибки на бит после декодирования по Хэммингу
Pe_bit_after_soft = zeros(1, length(SNRdB));		%Вероятность ошибки на бит после "мягкого декодера"
Ped = zeros(1,length(SNRdB));        				%Вероятность ошибки декодирования  
Ped_soft = zeros(1, length(SNRdB));					%Вероятность ошибки декодирования после "мягкого декодера"
T = zeros(1,length(SNRdB));  						%Пропускная способность канала

N = 100;
    
for i=1 : length(SNR)
    disp(i);
    sigma = sqrt(1/(2*SNR(i)));
    nTests = 0; nSent = 0;
    nErrDecode = 0; nErrDecode_soft = 0;
    nErrBits = 0; nErrBits_H = 0; nErrBits_H_soft = 0;
    messages_to_send = 0;
    while messages_to_send < N
        messages_to_send = messages_to_send + 1;
        ind = randi([1 2^k], 1, 1);												%Рандомим индекс
        c = reshape(mes_modBPSK(ind, :, :), 3, 12);								%Получаем сообщение по индеску
        while 1
            nTests = nTests + 1;
            r_noise = c + sigma*randn(3, 12);									%АБГШ
            unBPSK = r_noise < 0;												%Демодулятор BPSK
            nErrBits = nErrBits + sum(sum(xor(c<0, unBPSK)));
            plus_zeros = [unBPSK(:, :) [0 0 0; 0 0 0; 0 0 0]];					%Собираем 3 части в одно сообщение и добавляем 0
            corrected = zeros(3, 12);
            for part=1:3
                corr = plus_zeros(part, :); 									%Берем отдельную часть 
                S = bi2de(mod(corr*H',2)); 										%Получаем синдром отдельной части и производим синдромное декодирование
                switch S														%Для каждого синдрома посчитан вектор ошибки, который исправляет наиболее вероятную ошибку
                    case 1
                        corr = gfadd(corr,[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
                    case 2
                        corr = gfadd(corr,[0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);
                    case 3
                        corr = gfadd(corr,[0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]);
                    case 4
                        corr = gfadd(corr,[0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]);
                    case 5
                        corr = gfadd(corr,[0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]);
                    case 6
                        corr = gfadd(corr,[0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]);
                    case 7
                        corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]);
                    case 8
                        corr = gfadd(corr,[0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]);
                    case 9
                        corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]);
                    case 10
                        corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]);
                    case 11
                        corr = gfadd(corr,[0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]);
                    case 12
                        corr = gfadd(corr,[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]);
                    case 13
                        corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);
                    case 14
                        corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]);
                    case 15
                        corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]);
                end
                corrected(part, :) = corr(1, 1:end-3);		
            end
            nErrBits_H = nErrBits_H + sum(sum(xor(corrected, c<0)));
            concatenated = [corrected(1, 5:end) ...
                            corrected(2, 5:end) ...
                            corrected(3, 5:end)];								%Собираем все в одно сообщение, удаляя нулевые биты в начале
                        
                        
            [soft_corr, v_soft] = soft(r_noise, h, all_codewords, c);
            nErrBits_H_soft = nErrBits_H_soft + v_soft;
                        
                        
            [~, S] = gfdeconv(concatenated, g); 
            v = sum(xor(codes(ind, :), concatenated));
            
            [~, S_soft] = gfdeconv(soft_corr, g); 
            v_soft = sum(xor(codes(ind, :), soft_corr));
                       
            
            if (bi2de(S_soft) == 0) && (v_soft > 0)
                nErrDecode_soft = nErrDecode_soft + 1;
            end
            
            if (bi2de(S) == 0)
                if (v > 0)
                    nErrDecode = nErrDecode + 1;
                end
                nSent = nSent + 1;
                break;
            end
        end
    end
    Ped(i) = nErrDecode/nTests;
    Ped_soft(i) = nErrDecode_soft/nTests;
    Pe_bit(i) = nErrBits/(nTests*36);						%Делим на 36, т.к передавали 
    Pe_bit_after_soft(i) = nErrBits_H_soft / (nTests*36);
    Pe_bit_after_Hamming(i) = nErrBits_H/(nTests*36);
    T(i) = k * nSent / (36 * nTests);
end

nfig = 1;
figure(nfig)
semilogy(SNRdB, Ped, 'ko', ...
         SNRdB, Ped_soft, 'ro', ...
         SNRdB, Ped_theor_up, 'r.-', ...
         SNRdB, Ped_theor, 'b.-')
legend('Ped', 'Ped_soft', 'Ped theor up', 'Ped theor');


nfig = nfig + 1;
figure(nfig)
semilogy(SNRdB, Pe_bit, 'ko', ...
         SNRdB, Pe_bit_after_Hamming, 'b.--', ...
         SNRdB, Pe_bit_after_soft, 'r.--', ...
         SNRdB, Pe_bit_theor, 'm.-');
legend('Pe bit', 'Pe bit after Hamming', ...
       'Pe_bit_after_soft', 'Pe bit theor');
xlim([-10 0]);

nfig = nfig + 1;
figure(nfig)
plot(SNRdB, T, 'b.-');
ylabel('T, Пропускная способность');
xlabel('E/N_0, dB')


function [concatenated, v] = soft(AWGN, h, all_codewords, c)
    min_d = [100 100 100];
    corrected = all_codewords(1, :);
    cor = zeros(3, 12);
    for j=1:256
        for AWGN_part=1:3
            min = sqrt(sum((AWGN(AWGN_part, :) - h(j, :)).^2));
            if min < min_d(AWGN_part)
                min_d(AWGN_part) = min;
                corrected(AWGN_part, :) = all_codewords(j, :);
                cor(AWGN_part, :) = h(j, :);
            end
        end
    end
    v = sum(sum(xor(cor<0, c<0)));
    concatenated = [corrected(1, :) ...
                    corrected(2, :) ...
                    corrected(3, :)];
end

