% Alamouti 2x2 OFDM with LS Channel Estimation and MMSE Combining
% Flat Rayleigh fading, perfect synchronization

clear; clc;

% PARAMETERS
Nsub   = 64;
cpLen  = 16;
M      = 16;
bps    = log2(M);
SNRvec = 0:2:30;
numIter = 200;

rng(0);

BER_A22 = zeros(length(SNRvec),1);

% SNR LOOP
for snrIdx = 1:length(SNRvec)

    SNR = SNRvec(snrIdx);
    noiseVar = 10^(-SNR/10);

    bitErr = 0;
    bitTot = 0;

    for iter = 1:numIter

        %% DATA GENERATION
        bits = randi([0 1], 2*Nsub*bps, 1);
        bits_mat = reshape(bits, bps, []).';
        sym_idx = bi2de(bits_mat,'left-msb');
        sym = qammod(sym_idx, M, 'UnitAveragePower', true);

        s1 = sym(1:2:end);
        s2 = sym(2:2:end);

        
        pilot = qammod(zeros(Nsub,1), M, 'UnitAveragePower', true);

        % Pilot time slots
        X1_p1 = pilot;          X2_p1 = zeros(Nsub,1);
        X1_p2 = zeros(Nsub,1);  X2_p2 = pilot;

        
        X1_t1 = s1;           X2_t1 = s2;
        X1_t2 = -conj(s2);    X2_t2 = conj(s1);

        
        tx1_fd = [X1_p1 X1_p2 X1_t1 X1_t2];
        tx2_fd = [X2_p1 X2_p2 X2_t1 X2_t2];

        tx1 = ifft(tx1_fd, Nsub, 1);
        tx2 = ifft(tx2_fd, Nsub, 1);

        tx1 = [tx1(end-cpLen+1:end,:); tx1]; tx1 = tx1(:);
        tx2 = [tx2(end-cpLen+1:end,:); tx2]; tx2 = tx2(:);

        
        h11 = (randn + 1j*randn)/sqrt(2);
        h12 = (randn + 1j*randn)/sqrt(2);
        h21 = (randn + 1j*randn)/sqrt(2);
        h22 = (randn + 1j*randn)/sqrt(2);

        rx1 = awgn(h11*tx1 + h12*tx2, SNR);
        rx2 = awgn(h21*tx1 + h22*tx2, SNR);

        
        rx1_cp = reshape(rx1, Nsub+cpLen, []);
        rx2_cp = reshape(rx2, Nsub+cpLen, []);

        rx1_fd = fft(rx1_cp(cpLen+1:end,:), Nsub, 1);
        rx2_fd = fft(rx2_cp(cpLen+1:end,:), Nsub, 1);

        
        H11 = rx1_fd(:,1) ./ pilot;
        H12 = rx1_fd(:,2) ./ pilot;
        H21 = rx2_fd(:,1) ./ pilot;
        H22 = rx2_fd(:,2) ./ pilot;

        
        y11 = rx1_fd(:,3);  y12 = rx1_fd(:,4);
        y21 = rx2_fd(:,3);  y22 = rx2_fd(:,4);

        % Per-Rx Alamouti combining
        s1_r1 = conj(H11).*y11 + H12.*conj(y12);
        s2_r1 = conj(H12).*y11 - H11.*conj(y12);

        s1_r2 = conj(H21).*y21 + H22.*conj(y22);
        s2_r2 = conj(H22).*y21 - H21.*conj(y22);

        % Sum across receive antennas
        s1_hat = s1_r1 + s1_r2;
        s2_hat = s2_r1 + s2_r2;

        % MMSE normalization
        den = abs(H11).^2 + abs(H12).^2 + ...
              abs(H21).^2 + abs(H22).^2 + noiseVar;

        s1_hat = s1_hat ./ den;
        s2_hat = s2_hat ./ den;

        %% -------- DEMOD --------
        rx_sym = zeros(2*Nsub,1);
        rx_sym(1:2:end) = s1_hat;
        rx_sym(2:2:end) = s2_hat;

        rx_bits = reshape( ...
            de2bi(qamdemod(rx_sym, M, 'UnitAveragePower', true), ...
            bps, 'left-msb').', [], 1);

        bitErr = bitErr + sum(bits ~= rx_bits);
        bitTot = bitTot + length(bits);

    end

    BER_A22(snrIdx) = bitErr / bitTot;
    fprintf('SNR=%2d dB | Alamouti 2x2 MMSE BER = %.3e\n', ...
        SNR, BER_A22(snrIdx));
end

% ---------------- PLOT ----------------
figure;
semilogy(SNRvec, BER_A22, '-o','LineWidth',2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR: Alamouti 2x2 OFDM (LS + MMSE)');
