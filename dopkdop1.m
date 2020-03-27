clc
clear all
close all
nfig = 1;

k = 10;
g = [1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1];

gsize = length(g) - 1;
d_min = zeros(1, k);
x = zeros(1, k);

for j = 1:k
    n = j + gsize;
    codes = zeros(2^j, n);
    R = 2^gsize;
    for m = 0 : 2^j - 1
        mx = de2bi(m*R, n);
        mx = mx(end:-1:1);

        [~, r] = deconv(mx, g);
        cx = mod(r, 2);
        codes(m + 1, :) = xor(mx, cx);
    end

    d_min_cur = min(sum(codes(2:end, :), 2));
    d_min(j) = d_min_cur;
    x(j) = j;
end

figure()
plot(x, d_min);
axis([0 11 4 12]);
xlabel('word size');
ylabel('dmin');