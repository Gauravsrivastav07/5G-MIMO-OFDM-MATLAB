% SISO OFDM system simulation with AWGN channel

clear; clc;

Nsub = 64;
cpLen = 16;
M = 16;
SNR = 10;
bits_per_sym = log2(M);

%TRANSMITTER
bits = randi([0 1], Nsub*bits_per_sym, 1);
bits_mat = reshape(bits, bits_per_sym, []).';
sym_idx = bi2de(bits_mat, 'left-msb');
ofdm_fd = qammod(sym_idx, M, 'UnitAveragePower', true);

tx_time = ifft(ofdm_fd, Nsub);
tx_cp = [tx_time(end-cpLen+1:end); tx_time];
tx_serial = tx_cp(:);

% CHANNEL
rx_serial = awgn(tx_serial, SNR);

% RECEIVER
rx_cp_mat = reshape(rx_serial, Nsub+cpLen, []);
rx_no_cp = rx_cp_mat(cpLen+1:end, :);
rx_fd = fft(rx_no_cp, Nsub);
rx_sym = rx_fd(:);

rx_idx = qamdemod(rx_sym, M, 'UnitAveragePower', true);
rx_bits_mat = de2bi(rx_idx, bits_per_sym, 'left-msb');
rx_bits = reshape(rx_bits_mat.', [], 1);

BER = sum(bits ~= rx_bits) / length(bits)
