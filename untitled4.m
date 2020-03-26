clc
clear all
close all
nfig = 1;

k = 4;
N = 200000;
g = [1, 1, 1, 0, 1]; 
%g = [1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1];

gsize = length(g) - 1; 	%Получаем степень многочлена
n = k + gsize; 			%размер кодового слова 
codes = zeros(2^k, n); 	%кодовая книга
R = 2^gsize; 
for m = 0 : 2^k - 1
    mx = de2bi(m*R, n);	%Умножаем многочлен на x в степени нашего многочлена
    mx = mx(end:-1:1);	%Разворачиваем его
    
    [~, r] = deconv(mx, g);	%Получаем остаток от деления нашего многочлена на порождающий многочлен
    cx = mod(r, 2);			%Приводим остаток по mod 2
    codes(m + 1, :) = xor(mx, cx);	%Получаем кодовое слово
end

SNRdB = -20 : -10;
Peb = zeros(1, length(SNRdB));	%Массив для вероятности ошибки на бит
Ped = zeros(1, length(SNRdB));	%Массив для вероятности неправильного решения декодера
SNRtheor = 10.^(SNRdB/10);
Peb_theor = qfunc(sqrt(2*SNRtheor));	%Теоретическое значение вероятности ошибки на бит

for i = 1:length(SNRdB)
    SNR = 10.^(SNRdB(i)/10);
    sigma = sqrt(1/ (2*SNR));	
    errorBit = 0;	%Количество ошибок на бит для конкретного эксперимента
    errorDecode = 0;%Количество ошибок неправильного решения декодера для конкретного эксперимента
    for j = 1 : N
        ind = randi([1 2^k], 1, 1);	%Случайный индекс
        mx = codes(ind, :);			%Получаем кодовое слово по индексу
        r = mx.*-2 + 1;				%Модулируем по BPSK
        
        r_noise = r + sigma * randn(1, n);	%Добавляем шум
        
        mx_dec = r_noise < 0;		%Демодулируем
        
        e = xor(mx_dec, mx);			%Получаем вектор, в котором показывается сколько позиций не совпало(0 - совпало, 1 - не совпало)
        errorBit = errorBit + nnz(e);	%Увеличиваем количество ошибок на бит
        [~, r] = deconv(mx_dec, g);	
        cx = mod(r, 2);						%Получаем синдром
        if (nnz(e) > 0) & (sum(cx) == 0)
            errorDecode = errorDecode + 1;	%Увеличиваем количество ошибок неправильного решения
        end
    end
    Peb(1, i) = errorBit / N / n;	%Фиксируем результаты эксперимента
    Ped(1, i) = errorDecode / N;
end
        
Ped_as = (ones(1, length(SNRdB)) ./ 2).^ gsize;	%Асимптотическая верхняя граница вероятности ошибки декодирования


A = zeros(1, length(codes(1, :)) + 1);	%Создаем книгу весов
for i = 1: length(codes(:,1))
    w = sum(codes(i, :));
    A(1, w + 1) = A(1, w + 1) + 1;		
end

d_min = min(sum(codes(2:end, :), 2));	%Получаем минимальное расстояние кода
        
Ped_pr = zeros(1, length(SNRdB));
for i = 1 : length(SNRdB)
    for j = d_min : n
        Ped_pr(1, i) = Ped_pr(1, i) + A(j + 1) * Peb_theor(i)^j * (1 - Peb_theor(i))^(n - j);	%Получаем точное значение вероятности ошибки декодирования
    end
end

Ped_as2 = (2^k - 1).*Peb_theor.^d_min;	%Более точная верхняя граница вероятности ошибки декодирования

figure(nfig)
axis('square');
semilogy(SNRdB, Peb, 'rx', SNRdB, Peb_theor, 'g');
xlabel('SNRdB');
ylabel('PeBIT');
legend({'Практическое значение вероятности ошибки на бит', 'Теоретическое значение вероятности ошибки на бит'}, 'Location', 'SouthWest');

nfig = nfig + 1;
figure(nfig)
axis('square');
semilogy(SNRdB, Ped, 'bo', SNRdB, Ped_as, 'r', SNRdB, Ped_pr, 'm', SNRdB, Ped_as2, 'g');
xlabel('SNRdB');
ylabel('PeDECODE');

legend({'Практическое значение вероятности ошибки декодированиея',...
    'Асимптотическая верхняя граница вероятности ошибки декодирования', ...
    'Точное значение вероятности ошибки декодирования',...
    'Более точная верхняя граница вероятности ошибки декодирования'}, 'Location', 'SouthWest');

        
        
        
        
        
        
        
        