close all;
clear all;
clc;

%% ============================================================
% PERFECT CSI CHECK: OTFS and OFDM
%
% Purpose:
%   Compare G-matrix channel output and TDL/FIR channel output
%   under the same Perfect CSI receiver pipeline.
%
% Checks:
%   1. Clean channel error:
%        norm(r_TDL_clean - G_true*s) / norm(G_true*s)
%
%   2. Same-noise received error:
%        norm(r_TDL - r_G) / norm(r_G)
%
%   3. Equalized symbol error:
%        norm(estdata_TDL - estdata_G) / norm(estdata_G)
%
%   4. Decoded bit mismatch:
%        decoded_G vs decoded_TDL
%
% Runs:
%   Modulation      = OTFS, OFDM
%   Variant         = CP
%   CSI             = Perfect only
%   Detection       = LMMSE only
%   Channel coding  = None, Conv, Polar
%
% Main stabilizing changes:
%   1. polar_nMax_encode = polar_nMax_decode
%   2. exact LLR -> approxllr
%   3. LLR clipping before Conv / Polar decoding
% ============================================================

%% =========================
% Parameters
% =========================
N = 64;
M = 4;
MN = N*M;

M_mod  = 4;
M_bits = log2(M_mod);

eng_sqrt = (M_mod==2) + (M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

BandWidth = 2e6;
TdMax = 0.5e-6;
length_cp = ceil(TdMax*BandWidth);

data_grid = ones(M,N);

N_syms_perfram = sum(sum(data_grid));
N_bits_perfram = N_syms_perfram * M_bits;

car_fre = 18e9;
delta_f = BandWidth/M;
T = 1/delta_f;

payLoad_UncodedBits = 256;
codeRate = 0.5;

SNR_dB = 2:2:30;
SNR = 10.^(SNR_dB/10);

sigma_2_orig = (abs(eng_sqrt)^2)./SNR;

N_fram = 1000;
seedGen = 1:50:5000;

ModulationList = ["OTFS","OFDM"];
ChannelCodingList = ["None","Conv","Polar"];

Variant = "CP";
Detection = "LMMSE";
ChannelEstimation = "Perfect";

CheckTolerance = 1e-10;

% Polar consistency
polar_nMax_encode = 9;
polar_nMax_decode = 9;

polar_crcLen = 11;
polar_listSize = 8;

% LLR clipping
LLR_MAX = 50;

plot_floor = 1e-8;

%% =========================
% DFT matrices
% =========================
Fn = dftmtx(N);
Fn = Fn ./ norm(Fn);

Fm = dftmtx(M);
Fm = Fm ./ norm(Fm);

%% =========================
% Storage
% =========================
curveCount = 0;
curveLabels = {};

BER_G_curves = {};
BER_TDL_curves = {};
BER_Diff_curves = {};

FER_G_curves = {};
FER_TDL_curves = {};
FER_Diff_curves = {};

CleanErr_curves = {};
RxErr_curves = {};
SymErr_curves = {};
LLRErr_curves = {};

%% =========================
% Main sweep
% =========================
for Modulation = ModulationList

    for ChannelCoding = ChannelCodingList

        curveCount = curveCount + 1;
        curveLabels{curveCount} = sprintf('%s | %s', char(Modulation), char(ChannelCoding));

        errBits_G = zeros(1,length(SNR_dB));
        errBits_TDL = zeros(1,length(SNR_dB));
        diffBits_GT = zeros(1,length(SNR_dB));

        errFrames_G = zeros(1,length(SNR_dB));
        errFrames_TDL = zeros(1,length(SNR_dB));
        diffFrames_GT = zeros(1,length(SNR_dB));

        BER_G = zeros(1,length(SNR_dB));
        BER_TDL = zeros(1,length(SNR_dB));
        BER_Diff = zeros(1,length(SNR_dB));

        FER_G = zeros(1,length(SNR_dB));
        FER_TDL = zeros(1,length(SNR_dB));
        FER_Diff = zeros(1,length(SNR_dB));

        CleanErr = zeros(1,length(SNR_dB));
        RxErr = zeros(1,length(SNR_dB));
        SymErr = zeros(1,length(SNR_dB));
        LLRErr = zeros(1,length(SNR_dB));

        fprintf('\n============================================================\n');
        fprintf('Running: Modulation = %s | Coding = %s\n', Modulation, ChannelCoding);
        fprintf('============================================================\n');

        for iesn0 = 1:length(SNR_dB)

            noiseVar = sigma_2_orig(iesn0);
            frameCounter = 0;

            cleanErrSum = 0;
            rxErrSum = 0;
            symErrSum = 0;
            llrErrSum = 0;
            llrCount = 0;

            for a = 1:length(seedGen)

                seed = seedGen(a);

                for ifram = 1:N_fram

                    frameCounter = frameCounter + 1;

                    baseSeed  = 100000000 + 1000*a + ifram;
                    noiseSeed = 200000000 + 1000000*iesn0 + 1000*a + ifram;

                    rng(baseSeed);

                    %% ============================================
                    % TX frame generation
                    % ============================================
                    [trans_info_bit, X, meta] = make_tx_frame( ...
                        ChannelCoding, ...
                        N, M, M_mod, M_bits, data_grid, ...
                        N_syms_perfram, N_bits_perfram, ...
                        payLoad_UncodedBits, codeRate, ...
                        polar_nMax_encode, polar_crcLen);

                    %% ============================================
                    % Modulation: OTFS or OFDM
                    % ============================================
                    switch Modulation

                        case "OTFS"
                            X_tilda = X * Fn';

                        case "OFDM"
                            X_tilda = Fm' * X;

                        otherwise
                            error('Unknown modulation.');
                    end

                    s = reshape(X_tilda, MN, 1);

                    %% ============================================
                    % True channel
                    % ============================================
                    max_speed = 12000;

                    [chan_coef, delay_taps, Doppler_taps, taps] = ...
                        Generate_delay_Doppler_channel_parameters_local( ...
                            N, M, car_fre, delta_f, T, max_speed, seed, 1);

                    L_set = unique(delay_taps);
                    Lmax = max(delay_taps);

                    assert(length_cp >= Lmax, ...
                        'CP length must be >= maximum delay tap.');

                    [G_true, gs_true] = build_cp_G_true( ...
                        N, M, length_cp, chan_coef, delay_taps, Doppler_taps);

                    %% ============================================
                    % G clean output
                    % ============================================
                    r_G_clean = G_true * s;

                    %% ============================================
                    % TDL clean output
                    % ============================================
                    r_TDL_clean = build_cp_TDL_output( ...
                        s, N, M, length_cp, delta_f, T, ...
                        chan_coef, delay_taps, Doppler_taps);

                    %% ============================================
                    % Clean channel comparison
                    % ============================================
                    cleanErr = norm(r_TDL_clean - r_G_clean) / max(norm(r_G_clean), eps);
                    cleanErrSum = cleanErrSum + cleanErr;

                    if cleanErr > CheckTolerance
                        warning('Clean TDL/G mismatch: %s | %s | SNR=%d | frame=%d | err=%.6e', ...
                            Modulation, ChannelCoding, SNR_dB(iesn0), frameCounter, cleanErr);
                    end

                    %% ============================================
                    % Same noise for G and TDL
                    % ============================================
                    rng(noiseSeed);

                    noise_common = sqrt(noiseVar/2) * ...
                        (randn(MN,1) + 1i*randn(MN,1));

                    r_G = r_G_clean + noise_common;
                    r_TDL = r_TDL_clean + noise_common;

                    rxErr = norm(r_TDL - r_G) / max(norm(r_G), eps);
                    rxErrSum = rxErrSum + rxErr;

                    %% ============================================
                    % Perfect CSI detection and decoding
                    % Same gs_true passed to both
                    % ============================================
                    [decoded_G, estdata_G, rxSoft_G] = decode_perfect_cp_lmmse( ...
                        r_G, ChannelCoding, meta, ...
                        N, M, M_mod, noiseVar, data_grid, ...
                        gs_true, L_set, length_cp, Modulation, ...
                        N_bits_perfram, payLoad_UncodedBits, ...
                        polar_nMax_decode, polar_listSize, polar_crcLen, LLR_MAX);

                    [decoded_TDL, estdata_TDL, rxSoft_TDL] = decode_perfect_cp_lmmse( ...
                        r_TDL, ChannelCoding, meta, ...
                        N, M, M_mod, noiseVar, data_grid, ...
                        gs_true, L_set, length_cp, Modulation, ...
                        N_bits_perfram, payLoad_UncodedBits, ...
                        polar_nMax_decode, polar_listSize, polar_crcLen, LLR_MAX);

                    %% ============================================
                    % Symbol and LLR comparison
                    % ============================================
                    symErr = norm(estdata_TDL - estdata_G) / max(norm(estdata_G), eps);
                    symErrSum = symErrSum + symErr;

                    if ~isempty(rxSoft_G) && ~isempty(rxSoft_TDL)
                        llrErr = norm(rxSoft_TDL - rxSoft_G) / max(norm(rxSoft_G), eps);
                        llrErrSum = llrErrSum + llrErr;
                        llrCount = llrCount + 1;
                    end

                    %% ============================================
                    % Error counts
                    % ============================================
                    errors_G = sum(decoded_G ~= trans_info_bit);
                    errors_TDL = sum(decoded_TDL ~= trans_info_bit);
                    diff_GT = sum(decoded_G ~= decoded_TDL);

                    errBits_G(iesn0) = errBits_G(iesn0) + errors_G;
                    errBits_TDL(iesn0) = errBits_TDL(iesn0) + errors_TDL;
                    diffBits_GT(iesn0) = diffBits_GT(iesn0) + diff_GT;

                    if errors_G > 0
                        errFrames_G(iesn0) = errFrames_G(iesn0) + 1;
                    end

                    if errors_TDL > 0
                        errFrames_TDL(iesn0) = errFrames_TDL(iesn0) + 1;
                    end

                    if diff_GT > 0
                        diffFrames_GT(iesn0) = diffFrames_GT(iesn0) + 1;
                    end

                end
            end

            infoLen = length(trans_info_bit);

            BER_G(iesn0) = errBits_G(iesn0) / (frameCounter * infoLen);
            BER_TDL(iesn0) = errBits_TDL(iesn0) / (frameCounter * infoLen);
            BER_Diff(iesn0) = diffBits_GT(iesn0) / (frameCounter * infoLen);

            FER_G(iesn0) = errFrames_G(iesn0) / frameCounter;
            FER_TDL(iesn0) = errFrames_TDL(iesn0) / frameCounter;
            FER_Diff(iesn0) = diffFrames_GT(iesn0) / frameCounter;

            CleanErr(iesn0) = cleanErrSum / frameCounter;
            RxErr(iesn0) = rxErrSum / frameCounter;
            SymErr(iesn0) = symErrSum / frameCounter;

            if llrCount > 0
                LLRErr(iesn0) = llrErrSum / llrCount;
            else
                LLRErr(iesn0) = 0;
            end

            clc;
            disp('####################################################################')
            fprintf('Perfect CSI G vs TDL checker\n');
            fprintf('Modulation = %s | Variant = %s | Detection = %s | Coding = %s\n', ...
                Modulation, Variant, Detection, ChannelCoding);
            fprintf('SNR = %d dB | Frames = %d\n', SNR_dB(iesn0), frameCounter);
            fprintf('BER_G      = %.6e\n', BER_G(iesn0));
            fprintf('BER_TDL    = %.6e\n', BER_TDL(iesn0));
            fprintf('BER_G_vs_TDL_decoded_difference = %.6e\n', BER_Diff(iesn0));
            fprintf('FER_G      = %.6e\n', FER_G(iesn0));
            fprintf('FER_TDL    = %.6e\n', FER_TDL(iesn0));
            fprintf('FER_G_vs_TDL_decoded_difference = %.6e\n', FER_Diff(iesn0));
            fprintf('Avg clean channel error = %.6e\n', CleanErr(iesn0));
            fprintf('Avg noisy received error = %.6e\n', RxErr(iesn0));
            fprintf('Avg equalized symbol error = %.6e\n', SymErr(iesn0));
            fprintf('Avg LLR error = %.6e\n', LLRErr(iesn0));

        end

        BER_G_curves{curveCount} = BER_G;
        BER_TDL_curves{curveCount} = BER_TDL;
        BER_Diff_curves{curveCount} = BER_Diff;

        FER_G_curves{curveCount} = FER_G;
        FER_TDL_curves{curveCount} = FER_TDL;
        FER_Diff_curves{curveCount} = FER_Diff;

        CleanErr_curves{curveCount} = CleanErr;
        RxErr_curves{curveCount} = RxErr;
        SymErr_curves{curveCount} = SymErr;
        LLRErr_curves{curveCount} = LLRErr;

    end
end

%% =========================
% Plot 1: BER G vs TDL
% =========================
figure;
hold on;

for i = 1:curveCount
    semilogy(SNR_dB, max(BER_G_curves{i},plot_floor), '-o', 'LineWidth', 1.5);
    semilogy(SNR_dB, max(BER_TDL_curves{i},plot_floor), '--s', 'LineWidth', 1.5);
end

grid on;
grid minor;
xlabel('SNR (dB)');
ylabel('BER');
title('Perfect CSI: OTFS/OFDM BER Comparison, G vs TDL');

leg = {};
for i = 1:curveCount
    leg{end+1} = sprintf('%s | G', curveLabels{i});
    leg{end+1} = sprintf('%s | TDL', curveLabels{i});
end
legend(leg, 'Interpreter', 'none', 'Location', 'southwest');
ylim([plot_floor 1]);

%% =========================
% Plot 2: Decoded mismatch BER
% =========================
figure;
hold on;

for i = 1:curveCount
    semilogy(SNR_dB, max(BER_Diff_curves{i},plot_floor), '-o', 'LineWidth', 1.5);
end

grid on;
grid minor;
xlabel('SNR (dB)');
ylabel('Decoded bit mismatch rate');
title('Perfect CSI: Decoded Output Difference, G vs TDL');
legend(curveLabels, 'Interpreter', 'none', 'Location', 'best');
ylim([plot_floor 1]);

%% =========================
% Plot 3: Clean channel difference
% =========================
figure;
hold on;

for i = 1:curveCount
    semilogy(SNR_dB, max(CleanErr_curves{i},1e-18), '-o', 'LineWidth', 1.5);
end

grid on;
grid minor;
xlabel('SNR (dB)');
ylabel('norm(r_{TDL,clean} - Gs) / norm(Gs)');
title('Perfect CSI: Clean Channel Output Difference');
legend(curveLabels, 'Interpreter', 'none', 'Location', 'best');

%% =========================
% Plot 4: Equalized symbol difference
% =========================
figure;
hold on;

for i = 1:curveCount
    semilogy(SNR_dB, max(SymErr_curves{i},1e-18), '-o', 'LineWidth', 1.5);
end

grid on;
grid minor;
xlabel('SNR (dB)');
ylabel('norm(estdata_{TDL} - estdata_G) / norm(estdata_G)');
title('Perfect CSI: Equalized Symbol Difference');
legend(curveLabels, 'Interpreter', 'none', 'Location', 'best');

%% =========================
% Plot 5: LLR difference
% =========================
figure;
hold on;

for i = 1:curveCount
    semilogy(SNR_dB, max(LLRErr_curves{i},1e-18), '-o', 'LineWidth', 1.5);
end

grid on;
grid minor;
xlabel('SNR (dB)');
ylabel('norm(LLR_{TDL} - LLR_G) / norm(LLR_G)');
title('Perfect CSI: LLR Difference');
legend(curveLabels, 'Interpreter', 'none', 'Location', 'best');

%% =========================
% Final summary
% =========================
fprintf('\n============================================================\n');
fprintf('FINAL PERFECT CSI OTFS/OFDM G vs TDL SUMMARY\n');
fprintf('============================================================\n');

for i = 1:curveCount
    fprintf('\n%s\n', curveLabels{i});
    fprintf('Mean clean channel error  = %.6e\n', mean(CleanErr_curves{i}));
    fprintf('Mean noisy received error = %.6e\n', mean(RxErr_curves{i}));
    fprintf('Mean symbol error         = %.6e\n', mean(SymErr_curves{i}));
    fprintf('Mean LLR error            = %.6e\n', mean(LLRErr_curves{i}));
    fprintf('Mean decoded mismatch BER = %.6e\n', mean(BER_Diff_curves{i}));
end

fprintf('============================================================\n');

%% ============================================================
% Local functions
% ============================================================

function [trans_info_bit, X, meta] = make_tx_frame( ...
    ChannelCoding, ...
    N, M, M_mod, M_bits, data_grid, ...
    N_syms_perfram, N_bits_perfram, ...
    payLoad_UncodedBits, codeRate, ...
    polar_nMax_encode, polar_crcLen)

    meta = struct();
    meta.K_crc = [];
    meta.E = [];
    meta.crcLen = [];
    meta.puncturePattern = [];

    outLen = 512;

    if strcmp(ChannelCoding, "Conv")

        trans_info_bit = randi([0,1], payLoad_UncodedBits, 1);

        padLen = (outLen * codeRate) - payLoad_UncodedBits;

        if padLen ~= 0
            padBits = zeros(padLen, 1);
            mssgBits = [trans_info_bit; padBits];
        else
            mssgBits = trans_info_bit;
        end

        if codeRate == 0.5
            puncturePattern = [1 1];
        elseif codeRate == 0.75
            puncturePattern = [1 1 1 0 0 1];
        else
            error('Unsupported codeRate.');
        end

        trellis = poly2trellis(7, [171 133]);
        encodedBits = convenc(mssgBits, trellis, puncturePattern);

        if length(encodedBits) >= outLen
            codedBits = encodedBits(1:outLen);
        else
            reps = ceil(outLen / length(encodedBits));
            codedBits = repmat(encodedBits, reps, 1);
            codedBits = codedBits(1:outLen);
        end

        mappedBits = reshape(codedBits, M_bits, []);
        data = qammod(mappedBits, M_mod, 'gray', 'InputType', 'bit');

        X = Generate_2D_data_grid_local(N, M, data, data_grid);

        meta.puncturePattern = puncturePattern;

    elseif strcmp(ChannelCoding, "Polar")

        trans_info_bit = randi([0,1], payLoad_UncodedBits, 1);

        msgCRC = nrCRCEncode(trans_info_bit, '11');
        K_crc = length(msgCRC);

        E = N_syms_perfram * M_bits;
        iIL = false;

        enc = nrPolarEncode(msgCRC, E, polar_nMax_encode, iIL);
        modIn = nrRateMatchPolar(enc, K_crc, E, true);

        mappedBits = reshape(modIn, M_bits, []);
        data = qammod(mappedBits, M_mod, 'gray', 'InputType', 'bit');

        X = Generate_2D_data_grid_local(N, M, data, data_grid);

        meta.K_crc = K_crc;
        meta.E = E;
        meta.crcLen = polar_crcLen;

    elseif strcmp(ChannelCoding, "None")

        trans_info_bit = randi([0,1], N_bits_perfram, 1);

        data = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), ...
            M_mod, 'gray', 'InputType', 'bit');

        X = Generate_2D_data_grid_local(N, M, data, data_grid);

    else
        error('Unknown ChannelCoding option.');
    end
end

function X = Generate_2D_data_grid_local(N, M, data, data_grid)

    X_array = zeros(1, N*M);

    data_array = reshape(data_grid, 1, N*M);
    [~, data_index] = find(data_array > 0);

    data = data(:).';

    X_array(data_index) = data;

    X = reshape(X_array, M, N);
end

function [chan_coef,delay_taps,Doppler_taps,taps] = ...
    Generate_delay_Doppler_channel_parameters_local(N,M,car_fre,delta_f,T,max_speed,seed,U)

    rng(seed);

    one_delay_tap = 1/(M*delta_f*U);
    one_doppler_tap = 1/(N*T);

    K_factor = 0;

    delays = [0 310]*10^(-9);

    taps = length(delays);

    delay_taps = round(delays/one_delay_tap);

    pdp = [0 -3.6];

    pow_prof = 10.^(pdp/10);
    pow_prof = pow_prof/sum(pow_prof);

    rayleigh_part = sqrt(pow_prof).*(sqrt(1/2) * ...
        (randn(1,taps)+1i*randn(1,taps)));

    if K_factor > 0
        theta = 2*pi*rand(1,taps);
        los_part = sqrt(pow_prof).*exp(1i*theta);

        chan_coef = sqrt(K_factor/(K_factor+1))*los_part + ...
                    sqrt(1/(K_factor+1))*rayleigh_part;
    else
        chan_coef = rayleigh_part;
    end

    max_UE_speed = max_speed*(1000/3600);
    Doppler_vel = (max_UE_speed*car_fre)/(299792458);
    max_Doppler_tap = Doppler_vel/one_doppler_tap;

    Doppler_taps = max_Doppler_tap*cos(2*pi*rand(1,taps));
end

function [G_true, gs_true] = build_cp_G_true( ...
    N, M, length_cp, chan_coef, delay_taps, Doppler_taps)

    MN = N*M;

    l_i_true = delay_taps(:).';
    g_i_true = chan_coef(:).';
    k_i_true = Doppler_taps(:).';

    taps = length(l_i_true);
    Lmax = max(l_i_true);

    z = exp(1i*2*pi/MN);

    gs_true = zeros(Lmax+1, N*(M+length_cp));

    for q = 0:N*(M+length_cp)-1
        for i = 1:taps
            gs_true(l_i_true(i)+1, q+1) = gs_true(l_i_true(i)+1, q+1) + ...
                g_i_true(i) * z^(k_i_true(i)*(q-l_i_true(i)));
        end
    end

    G_true = zeros(MN, MN);

    for n = 0:N-1
        for m = 0:M-1
            for ell = 0:Lmax
                G_true(m+n*M+1, n*M+mod(m-ell,M)+1) = ...
                    gs_true(ell+1, m+n*(M+length_cp)+1);
            end
        end
    end
end

function r_clean = build_cp_TDL_output( ...
    s, N, M, length_cp, delta_f, T, ...
    chan_coef, delay_taps, Doppler_taps)

    MN = N*M;

    l_i_true = delay_taps(:).';
    g_i_true = chan_coef(:).';
    k_i_true = Doppler_taps(:).';

    taps = length(l_i_true);
    Lmax = max(l_i_true);

    offset = length_cp;

    s_ext = zeros(N*(M+length_cp), 1);

    for n = 0:N-1
        block_s = s(n*M+1 : (n+1)*M);
        s_ext(n*(M+length_cp)+1 : (n+1)*(M+length_cp)) = ...
            [block_s(end-length_cp+1:end); block_s];
    end

    L_in = length(s_ext);

    Fs = M * delta_f;
    Doppler_resolution = 1 / (N * T);
    f_d_taps = k_i_true * Doppler_resolution;

    gs_f = zeros(Lmax + 1, L_in);

    for q = 0:(L_in - 1)

        q_eff = q - offset;

        for i = 1:taps

            ell = l_i_true(i);

            physical_phase = exp(1i * 2*pi * f_d_taps(i) * (q_eff - ell) / Fs);

            gs_f(ell + 1, q + 1) = gs_f(ell + 1, q + 1) + ...
                g_i_true(i) * physical_phase;

        end
    end

    r_fir = zeros(L_in, 1);

    for q = 0:(L_in - 1)
        for ell = 0:Lmax
            if q >= ell
                r_fir(q + 1) = r_fir(q + 1) + ...
                    gs_f(ell + 1, q + 1) * s_ext(q - ell + 1);
            end
        end
    end

    r_clean = zeros(MN, 1);

    for n = 0:N-1

        idx_start = n*(M+length_cp) + length_cp + 1;
        idx_end = (n+1)*(M+length_cp);

        r_clean(n*M+1 : (n+1)*M) = r_fir(idx_start : idx_end);

    end
end

function [decoded_Bits, estdata, rxSoft] = decode_perfect_cp_lmmse( ...
    r, ChannelCoding, meta, ...
    N, M, M_mod, noiseVar, data_grid, ...
    gs_true, L_set, length_cp, Modulation, ...
    N_bits_perfram, payLoad_UncodedBits, ...
    polar_nMax_decode, polar_listSize, polar_crcLen, LLR_MAX)

    rxSoft = [];

    [estBits, estdata] = Block_LMMSE_detector_local( ...
        N, M, M_mod, noiseVar, data_grid, r, ...
        gs_true, L_set, length_cp, 1, Modulation);

    if strcmp(ChannelCoding, "None")

        decoded_Bits = estBits;

    elseif strcmp(ChannelCoding, "Conv")

        trellis = poly2trellis(7, [171 133]);

        rxSoft = make_clipped_llr(estdata, M_mod, noiseVar, N_bits_perfram, LLR_MAX);

        traceback = 35;

        decoded_Bits = vitdec(rxSoft, trellis, traceback, ...
            'trunc', 'unquant', meta.puncturePattern);

        decoded_Bits = decoded_Bits(1:payLoad_UncodedBits);

    elseif strcmp(ChannelCoding, "Polar")

        rxSoft = make_clipped_llr(estdata, M_mod, noiseVar, N_bits_perfram, LLR_MAX);

        iIL = false;
        iBIL = true;

        recPolar = nrRateRecoverPolar(rxSoft, meta.K_crc, meta.E, iBIL);

        decoded_crc = nrPolarDecode( ...
            recPolar, ...
            meta.K_crc, ...
            meta.E, ...
            polar_listSize, ...
            polar_nMax_decode, ...
            iIL, ...
            polar_crcLen);

        decoded_Bits = decoded_crc(1:payLoad_UncodedBits);

    else
        error('Unknown ChannelCoding option.');
    end
end

function rxSoft = make_clipped_llr(estdata, M_mod, noiseVar, N_bits_perfram, LLR_MAX)

    rxSoft = reshape(qamdemod(estdata, M_mod, 'gray', ...
        'OutputType', 'approxllr', ...
        'NoiseVariance', noiseVar), ...
        N_bits_perfram, 1);

    rxSoft(isnan(rxSoft)) = 0;
    rxSoft(isinf(rxSoft) & rxSoft > 0) = LLR_MAX;
    rxSoft(isinf(rxSoft) & rxSoft < 0) = -LLR_MAX;

    rxSoft(rxSoft >  LLR_MAX) =  LLR_MAX;
    rxSoft(rxSoft < -LLR_MAX) = -LLR_MAX;
end

function [est_bits,x_data] = Block_LMMSE_detector_local( ...
    N,M,M_mod,noise_var,data_grid,r,gs,L_set,length_cp,U,Modulation)

    Fn = dftmtx(N);
    Fn = Fn./norm(Fn);

    Fm = dftmtx(M);
    Fm = Fm./norm(Fm);

    N_syms_perfram = sum(sum((data_grid > 0)));

    data_array = reshape(data_grid,1,N*M);
    [~,data_index] = find(data_array > 0);

    M_bits = log2(M_mod);
    N_bits_perfram = N_syms_perfram*M_bits;

    sn_block_est = zeros(M,N);

    for n = 1:N

        Gn = zeros(M*U,M*U);

        for m = 1:M*U
            for l = L_set
                Gn(m,mod(m-l-1,M*U)+1) = ...
                    gs(l+1,(m-1)+ (n-1)*(M+length_cp)*U+1);
            end
        end

        rn = r((n-1)*M*U+1:n*M*U);

        Rn = Gn' * Gn;

        sn_hat = (Rn + noise_var.*eye(M*U)) \ (Gn' * rn);

        if U == 1
            sn_block_est(:,n) = sn_hat;
        else
            sn_block_est(:,n) = resample(sn_hat,1,U);
        end
    end

    X_tilda_est = sn_block_est;

    if strcmp(Modulation,"OTFS")
        X_est = X_tilda_est * Fn;
    elseif strcmp(Modulation,"OFDM")
        X_est = Fm * X_tilda_est;
    else
        error('Unknown Modulation.');
    end

    x_est = reshape(X_est,1,N*M);
    x_data = x_est(data_index);

    est_bits = reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'), ...
        N_bits_perfram,1);
end