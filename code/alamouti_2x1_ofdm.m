% Alamouti 2x1 OFDM with LS channel estimation (flat Rayleigh fading)

clc; clear;

Nsub = 64;
cpLen = 16;
M = 16;
bps = log2(M);
SNR = 10;
numIter = 500;

err = 0; totBits = 0;

for iter = 1:numIter

    bits = randi([0 1], 2*Nsub*bps, 1);
    bits_mat = reshape(bits, bps, []).';
    bits_idx = bi2de(bits_mat, 'left-msb');
    sym = qammod(bits_idx, M, 'UnitAveragePower', true);

    s1 = sym(1:2:end);
    s2 = sym(2:2:end);

    pilot = qammod(zeros(Nsub,1), M, 'UnitAveragePower', true);

    x1_p1 = pilot;      x2_p1 = zeros(Nsub,1);
    x1_p2 = zeros(Nsub,1);  x2_p2 = pilot;

    x1_t1 = s1;         x2_t1 = s2;
    x1_t2 = -conj(s2);  x2_t2 = conj(s1);

    tx1_fd = [x1_p1 x1_p2 x1_t1 x1_t2];
    tx2_fd = [x2_p1 x2_p2 x2_t1 x2_t2];

    tx1_time = ifft(tx1_fd, Nsub, 1);
    tx2_time = ifft(tx2_fd, Nsub, 1);

    tx1_cp = [tx1_time(end-cpLen+1:end,:); tx1_time];
    tx2_cp = [tx2_time(end-cpLen+1:end,:); tx2_time];

    tx1 = tx1_cp(:);
    tx2 = tx2_cp(:);

    % Flat Rayleigh fading (constant over subcarriers)
    h1 = (randn + 1j*randn)/sqrt(2);
    h2 = (randn + 1j*randn)/sqrt(2);

    rx = awgn(h1*tx1 + h2*tx2, SNR);

    rx_cp = reshape(rx, Nsub+cpLen, []);
    rx_no_cp = rx_cp(cpLen+1:end,:);
    rx_fd = fft(rx_no_cp, Nsub, 1);

    yp1 = rx_fd(:,1);
    yp2 = rx_fd(:,2);

    h1_est = yp1 ./ pilot;
    h2_est = yp2 ./ pilot;

    y1 = rx_fd(:,3);
    y2 = rx_fd(:,4);

    s1_hat = conj(h1_est).*y1 + h2_est.*conj(y2);
    s2_hat = conj(h2_est).*y1 - h1_est.*conj(y2);

    den = abs(h1_est).^2 + abs(h2_est).^2;
    s1_hat = s1_hat ./ den;
    s2_hat = s2_hat ./ den;

    rx_sym = zeros(2*Nsub,1);
    rx_sym(1:2:end) = s1_hat;
    rx_sym(2:2:end) = s2_hat;

    rx_idx = qamdemod(rx_sym, M, 'UnitAveragePower', true);
    rx_bits = reshape(de2bi(rx_idx, bps, 'left-msb').', [], 1);

    err = err + sum(rx_bits ~= bits);
    totBits = totBits + length(bits);
end

BER = err/totBits;
fprintf('Alamouti 2x1 OFDM (LS Estimation) | BER = %.3e\n', BER);
