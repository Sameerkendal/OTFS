close all
clear all
clc

%% ============================================================
% Merged RCP_OTFS.m
%
% Supported ChannelEstimation modes:
%   "Perfect"
%   "Estimator"
%   "Noisy"
%   "SP"          % superimposed pilots, supports "None", "Polar", "Conv"
%
% Supported ChannelCoding modes:
%   "Conv"
%   "Polar"
%   "None"
%
% Supported Detection modes:
%   "LMMSE"
%   "MRC"
%   "SP_LMMSE"    % used only when ChannelEstimation = "SP"
%
% Notes:
% 1) Original RCP paths ("Perfect", "Estimator", "Noisy") are kept.
% 2) SP mode is enabled for "None", "Polar", and "Conv".
% 3) When ChannelEstimation = "SP", the script automatically forces:
%       ChannelCodingList ⊂ {"Conv","Polar","None"}
%       DetectionList     = "SP_LMMSE"
% 4) SP + Polar / SP + Conv use approximate soft LLRs from post-equalized symbols.
% ============================================================

%% =========================
% Plot styling
% =========================
colors = {'r','g','c','m','y','k','b',[0.85 0.33 0.10],[0.49 0.18 0.56],[0.47 0.67 0.19]};
marker = {'o','+','*','s','d','^','x','p','v','>'};

%% =========================
% OTFS parameters
% =========================
N = 64;
M = 4;
MN = N*M;

M_mod  = 4;
M_bits = log2(M_mod);

eng_sqrt = (M_mod==2) + (M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

BandWidth = 2e6;
TdMax = 0.5e-6; 

data_grid = ones(M,N);

N_syms_perfram = sum(sum(data_grid));
N_bits_perfram = N_syms_perfram * M_bits;

car_fre = 18e9;
delta_f = BandWidth / M;
T = 1 / delta_f;

payLoad_UncodedBits = 384;
codeRate = 0.75;

SNR_dB = 0:5:20;
SNR    = 10.^(SNR_dB/10);

% Original RCP noise model (for Perfect / Estimator / Noisy)
sigma_2_orig = (abs(eng_sqrt)^2) ./ SNR;

% SP branch uses UnitAveragePower symbols, so use 1/SNR there
sigma_2_sp = 1 ./ SNR;

N_fram  = 10;
seedGen = 1:50:5000;

%% =========================
% User switches
% =========================
ChannelEstimation = ["Perfect","SP","Estimator"]          ; %"Perfect", "Estimator", "Noisy", "SP"
Modulation        = 0;             % 0 = OTFS, 1 = OFDM
ChannelCodingList = ["Conv","Polar","None"];
DetectionList     = ["LMMSE"];     % for non-SP: "LMMSE" or "MRC"

% SP parameters
Pd = 0.70;
Pp = 0.30;
assert(abs(Pd + Pp - 1) < 1e-12, 'Pd + Pp must equal 1.');

lambda_h     = 1e-4;
lambda_x     = 1e-3;
convTol_h    = 1e-5;
numOuterIter = 6;

% Automatic guard for SP mode
if strcmp(ChannelEstimation, "SP")
    allowedSPCoding = ["Conv","Polar","None"];
    keepMask = ismember(ChannelCodingList, allowedSPCoding);
    ChannelCodingList = ChannelCodingList(keepMask);

    if isempty(ChannelCodingList)
        ChannelCodingList = allowedSPCoding;
    end

    DetectionList = "SP_LMMSE";
    disp('SP mode selected: allowing only ChannelCodingList = ["Conv","Polar","None"] and DetectionList = "SP_LMMSE".');
end

%% =========================
% Normalized DFT matrices
% =========================
Fn = dftmtx(N);
Fn = Fn ./ norm(Fn);

Fm = dftmtx(M);
Fm = Fm ./ norm(Fm);

%% =========================
% Fixed known pilot pattern for SP mode
% =========================
rng(999);
pilotIdxFixed = randi([0, M_mod-1], N_syms_perfram, 1);
pilotSymFixed = qammod(pilotIdxFixed, M_mod, 'gray', ...
    'InputType', 'integer', 'UnitAveragePower', true);
Xpilot_unit_fixed = Generate_2D_data_grid(N, M, pilotSymFixed, data_grid);
Xpilot_fixed = sqrt(Pp) * Xpilot_unit_fixed;

%% =========================
% Storage for all curves
% =========================
curveCount  = 0;
curveLabels = {};
FER_curves  = {};
BER_curves  = {};
NMSE_curves = {};

%% =========================
% Main sweep
% =========================
for ChannelCoding = ChannelCodingList
    for Detection = DetectionList

        % Safety guards
        if strcmp(ChannelEstimation, "SP") && ~any(strcmp(ChannelCoding, ["None","Polar","Conv"]))
            continue;
        end
        if ~strcmp(ChannelEstimation, "SP") && strcmp(Detection, "SP_LMMSE")
            continue;
        end
        if strcmp(ChannelEstimation, "SP") && ~strcmp(Detection, "SP_LMMSE")
            continue;
        end

        curveCount = curveCount + 1;

        current_frame_number = zeros(1, length(SNR_dB));

        err_ber = zeros(1, length(SNR_dB));
        avg_ber = zeros(1, length(SNR_dB));

        err_fer = zeros(1, length(SNR_dB));
        avg_fer = zeros(1, length(SNR_dB));

        nmse_sum = zeros(1, length(SNR_dB));
        avg_nmse = zeros(1, length(SNR_dB));

        for iesn0 = 1:length(SNR_dB)

            frameCounter = 0;

            for a = 1:length(seedGen)
                seed = seedGen(a);

                for ifram = 1:N_fram

                    frameCounter = frameCounter + 1;
                    current_frame_number(iesn0) = ifram;

                    % reset frame-specific vars
                    K_crc = [];
                    crcLen = [];
                    E = [];
                    puncturePattern = [];
                    Xdata = [];
                    Xpilot = [];
                    s_pilot = [];

                    %% ============================================
                    % FRAME GENERATION / ENCODING
                    % ============================================
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

                        if(codeRate == 0.5)
                            puncturePattern = [1 1];
                        elseif(codeRate == 0.75)
                            puncturePattern = [1 1 1 0 0 1];
                        else
                            error('Unsupported codeRate for Conv branch.');
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

                        if strcmp(ChannelEstimation, "SP")
                            dataSym = qammod(mappedBits, M_mod, 'gray', ...
                                'InputType', 'bit', 'UnitAveragePower', true);

                            Xdata_unit = Generate_2D_data_grid(N, M, dataSym, data_grid);
                            Xdata = sqrt(Pd) * Xdata_unit;

                            Xpilot = Xpilot_fixed;
                            X = Xdata + Xpilot;
                        else
                            data = qammod(mappedBits, M_mod, 'gray', 'InputType', 'bit');
                            X = Generate_2D_data_grid(N, M, data, data_grid);
                        end

                    elseif strcmp(ChannelCoding, "Polar")

                        trans_info_bit = randi([0,1], payLoad_UncodedBits, 1);

                        mssgBits = trans_info_bit;
                        crcLen = 11;
                        msgCRC = nrCRCEncode(mssgBits, '11');
                        K_crc = length(msgCRC);
                        E = N_syms_perfram * M_bits;
                        nMax = 9;
                        iIL = false;

                        enc = nrPolarEncode(msgCRC, E, nMax, iIL);
                        modIn = nrRateMatchPolar(enc, length(msgCRC), E, true);

                        mappedBits = reshape(modIn, M_bits, []);

                        if strcmp(ChannelEstimation, "SP")
                            dataSym = qammod(mappedBits, M_mod, 'gray', ...
                                'InputType', 'bit', 'UnitAveragePower', true);

                            Xdata_unit = Generate_2D_data_grid(N, M, dataSym, data_grid);
                            Xdata = sqrt(Pd) * Xdata_unit;

                            Xpilot = Xpilot_fixed;
                            X = Xdata + Xpilot;
                        else
                            data = qammod(mappedBits, M_mod, 'gray', 'InputType', 'bit');
                            X = Generate_2D_data_grid(N, M, data, data_grid);
                        end

                    elseif strcmp(ChannelCoding, "None")

                        if strcmp(ChannelEstimation, "SP")
                            trans_info_bit = randi([0,1], N_bits_perfram, 1);

                            dataSym = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), ...
                                             M_mod, 'gray', 'InputType', 'bit', ...
                                             'UnitAveragePower', true);

                            Xdata_unit = Generate_2D_data_grid(N, M, dataSym, data_grid);
                            Xdata = sqrt(Pd) * Xdata_unit;

                            Xpilot = Xpilot_fixed;
                            X = Xdata + Xpilot;
                        else
                            trans_info_bit = randi([0,1], N_bits_perfram, 1);
                            data = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), ...
                                          M_mod, 'gray', 'InputType', 'bit');
                            X = Generate_2D_data_grid(N, M, data, data_grid);
                        end

                    else
                        error('Unknown ChannelCoding option.');
                    end

                    %% ============================================
                    % MODULATION
                    % ============================================
                    if strcmp(ChannelEstimation, "SP") && any(strcmp(ChannelCoding, ["None","Polar","Conv"]))

                        if(Modulation == 0)   % OTFS
                            Xdata_tilda  = Xdata * Fn';
                            Xpilot_tilda = Xpilot * Fn';
                            X_tilda      = X * Fn';
                        elseif(Modulation == 1) % OFDM
                            Xdata_tilda  = Fm' * Xdata;
                            Xpilot_tilda = Fm' * Xpilot;
                            X_tilda      = Fm' * X;
                        else
                            error('Unknown Modulation option.');
                        end

                        s_data  = reshape(Xdata_tilda,  MN, 1); 
                        s_pilot = reshape(Xpilot_tilda, MN, 1);
                        s       = reshape(X_tilda,      MN, 1);

                        noiseVar = sigma_2_sp(iesn0);

                    else

                        if(Modulation == 0)   % OTFS
                            X_tilda = X * Fn';
                        elseif(Modulation == 1) % OFDM
                            X_tilda = Fm' * X;
                        else
                            error('Unknown Modulation option.');
                        end

                        s = reshape(X_tilda, MN, 1);
                        noiseVar = sigma_2_orig(iesn0);
                    end

                    %% ============================================
                    % TRUE CHANNEL CONSTRUCTION
                    % ============================================
                    max_speed = 12000;

                    [chan_coef, delay_taps, Doppler_taps, taps] = ...
                        Generate_delay_Doppler_channel_parameters( ...
                            N, M, car_fre, delta_f, T, max_speed, seed, 1);

                    L_set = unique(delay_taps);
                    Lmax = max(delay_taps);

                    gs_true = zeros(Lmax+1, MN);
                    z = exp(1i*2*pi/MN);

                    l_i_true = delay_taps(:).';
                    g_i_true = chan_coef(:).';
                    k_i_true = Doppler_taps(:).';

                    for q = 0:MN-1
                        for i = 1:taps
                            gs_true(l_i_true(i)+1, q+1) = gs_true(l_i_true(i)+1, q+1) + ...
                                g_i_true(i) * z^(k_i_true(i) * (q - l_i_true(i)));
                        end
                    end

                    G_rcp_true = zeros(MN, MN);
                    for q = 0:MN-1
                        for ell = 0:Lmax
                            G_rcp_true(q+1, mod(q-ell, MN)+1) = gs_true(ell+1, q+1);
                        end
                    end

                    %% ============================================
                    % CHANNEL OUTPUT
                    % ============================================
                    r = G_rcp_true * s;
                    noise = sqrt(noiseVar/2) * (randn(size(r)) + 1i*randn(size(r)));
                    r = r + noise;

                    %% ============================================
                    % DEMODULATED OBSERVATION
                    % ============================================
                    Y_tilda = reshape(r, M, N);
                    if(Modulation == 0)
                        Y = Y_tilda * Fn;
                    elseif(Modulation == 1)
                        Y = Fm * Y_tilda;
                    else
                        error('Unknown Modulation option.');
                    end

                    %% ============================================
                    % CHANNEL ESTIMATION / CHANNEL MODEL TO USE
                    % ============================================
                    G_used = G_rcp_true;
                    gs_used = gs_true;
                    estBits = [];
                    estdata = [];
                    decoded_Bits = [];
                    softSymUnit = [];

                    switch ChannelEstimation

                        case "Perfect"

                            G_used = G_rcp_true;
                            gs_used = gs_true;

                        case "Estimator"

                            xp = 6;
                            pilot_grid = zeros(M, N);
                            pilot_grid(2, 32) = 1;
                            pilot_positions = find(pilot_grid);

                            N_pilot_syms = sum(sum(pilot_grid));
                            pilot_data = xp * ones(N_pilot_syms, 1);
                            X_pilot = Generate_2D_data_grid(N, M, pilot_data, pilot_grid);

                            X_pilot_tilda = X_pilot * Fn';
                            s_pilot_est = reshape(X_pilot_tilda, MN, 1);

                            r_pilot = G_rcp_true * s_pilot_est;
                            noise_pilot = sqrt(noiseVar/2) * ...
                                (randn(size(s_pilot_est)) + 1i*randn(size(s_pilot_est)));
                            r_pilot = r_pilot + noise_pilot;

                            Y_pilot_tilda = reshape(r_pilot, M, N);
                            Y_pilot = Y_pilot_tilda * Fn;

                            max_delay_taps = max(delay_taps);
                            max_doppler_taps = 31;

                            [pilot_m, pilot_n] = ind2sub([M, N], pilot_positions(1));

                            k_est_range = max(1, pilot_n - max_doppler_taps) : min(N, pilot_n + max_doppler_taps);
                            l_est_range = max(1, pilot_m) : min(M, pilot_m + max_delay_taps);

                            threshold = 5 * sqrt(noiseVar);

                            l_i_est = [];
                            k_i_est = [];
                            g_i_est = [];
                            count = 1;

                            for k = k_est_range
                                for l = l_est_range
                                    y_received = Y_pilot(l, k);

                                    if abs(y_received) > threshold
                                        k_rel = k - pilot_n;
                                        l_rel = l - pilot_m;
                                        phase_compensation = exp(-1j * 2 * pi * (1/(M*N)) * k_rel * (pilot_m-1));

                                        l_i_est(count) = l_rel; 
                                        k_i_est(count) = k_rel; 
                                        g_i_est(count) = (y_received / xp) * phase_compensation;
                                        count = count + 1;
                                    end
                                end
                            end

                            taps_est = length(l_i_est);
                            gs_est = zeros(Lmax+1, MN);

                            for q = 0:MN-1
                                for i = 1:taps_est
                                    gs_est(l_i_est(i)+1, q+1) = gs_est(l_i_est(i)+1, q+1) + ...
                                        g_i_est(i) * z^(k_i_est(i) * (q - l_i_est(i)));
                                end
                            end

                            G_rcp_est = zeros(MN, MN);
                            for q = 0:MN-1
                                for ell = 0:Lmax
                                    G_rcp_est(q+1, mod(q-ell, MN)+1) = gs_est(ell+1, q+1);
                                end
                            end

                            G_used = G_rcp_est;
                            gs_used = gs_est;

                        case "Noisy"

                            EstmationGain_dB = 10;
                            pilotSNR_dB = SNR_dB(iesn0) + EstmationGain_dB;
                            pilotSNR = 10^(pilotSNR_dB/10);

                            Ppilot = 1;
                            Npilots = 1;
                            sigma_e2 = (Ppilot / pilotSNR) / (abs(Ppilot)^2 * Npilots);

                            e_i = sqrt(sigma_e2/2) * (randn(size(g_i_true)) + 1i*randn(size(g_i_true)));
                            g_i_hat = g_i_true + e_i;

                            gs_hat = zeros(Lmax+1, MN);
                            for q = 0:MN-1
                                for i = 1:taps
                                    gs_hat(l_i_true(i)+1, q+1) = gs_hat(l_i_true(i)+1, q+1) + ...
                                        g_i_hat(i) * z^(k_i_true(i) * (q - l_i_true(i)));
                                end
                            end

                            G_rcp_hat = zeros(MN, MN);
                            for q = 0:MN-1
                                for ell = 0:Lmax
                                    G_rcp_hat(q+1, mod(q-ell, MN)+1) = gs_hat(ell+1, q+1);
                                end
                            end

                            G_used = G_rcp_hat;
                            gs_used = gs_hat;

                        case "SP"

                            if ~any(strcmp(ChannelCoding, ["None","Polar","Conv"]))
                                error('SP mode is supported only for ChannelCoding = "None", "Polar", or "Conv".');
                            end

                            Psp = length(l_i_true);

                            [A_pilot, Gbasis] = build_sp_regression_matrix_rcp( ...
                                N, M, l_i_true, k_i_true, s_pilot);

                            hHat = (A_pilot' * A_pilot + lambda_h * eye(Psp)) \ (A_pilot' * r);

                            G_hat = zeros(MN, MN);
                            for p = 1:Psp
                                G_hat = G_hat + hHat(p) * Gbasis{p};
                            end

                            [estBits, estdata, s_data_hat, softSymUnit] = sp_lmmse_detector_rcp( ...
                                N, M, M_mod, noiseVar, r, G_hat, s_pilot, Pd, Modulation, lambda_x);

                            for it = 1:numOuterIter
                                hPrev = hHat;

                                s_ref = s_pilot + s_data_hat;

                                [A_ref, Gbasis] = build_sp_regression_matrix_rcp( ...
                                    N, M, l_i_true, k_i_true, s_ref);

                                hHat = (A_ref' * A_ref + lambda_h * eye(Psp)) \ (A_ref' * r);

                                G_hat = zeros(MN, MN);
                                for p = 1:Psp
                                    G_hat = G_hat + hHat(p) * Gbasis{p};
                                end

                                [estBits, estdata, s_data_hat, softSymUnit] = sp_lmmse_detector_rcp( ...
                                    N, M, M_mod, noiseVar, r, G_hat, s_pilot, Pd, Modulation, lambda_x);

                                deltaH = norm(hHat - hPrev)^2 / max(norm(hPrev)^2, eps);
                                if deltaH < convTol_h
                                    break;
                                end
                            end

                            G_used = G_hat;
                            gs_used = [];

                        otherwise
                            error('Unknown ChannelEstimation option.');
                    end

                    %% ============================================
                    % CHANNEL NMSE
                    % ============================================
                    nmse_this_frame = norm(G_used - G_rcp_true, 'fro')^2 / ...
                                      max(norm(G_rcp_true, 'fro')^2, eps);
                    nmse_sum(iesn0) = nmse_sum(iesn0) + nmse_this_frame;
                    avg_nmse(iesn0) = nmse_sum(iesn0) / frameCounter;

                    %% ============================================
                    % DETECTION / DECODING
                    % ============================================
                    if ~strcmp(ChannelEstimation, "SP")

                        if strcmp(Detection, "MRC")

                            n_ite_MRC = 100;
                            omega = 1;
                            if(M_mod == 64)
                                omega = 0.25;
                            end
                            decision = 0;
                            init_estimate = 0;

                            [~, ~, K_ml] = Gen_DT_and_DD_channel_vectors(N, M, Lmax, gs_used);
                            H_tf = Generate_time_frequency_channel(N, M, G_used);

                            [estBits, ~, estdata] = MRC_delaydoppler( ...
                                N, M, M, M_mod, noiseVar, data_grid, r, H_tf, ...
                                L_set, omega, decision, init_estimate, n_ite_MRC, K_ml, Y);

                        elseif strcmp(Detection, "LMMSE")

                            [estBits, estdata] = LMMSE_detector( ...
                                N, M, M_mod, noiseVar, data_grid, r, G_used, Modulation);

                        else
                            error('Unsupported Detection option for non-SP mode.');
                        end

                        if strcmp(ChannelCoding, "Conv")

                            trellis = poly2trellis(7, [171 133]);
                            rxSoft = reshape(qamdemod(estdata, M_mod, 'gray', ...
                                'OutputType', 'llr', 'NoiseVariance', noiseVar), ...
                                N_bits_perfram, 1);

                            traceback = 35;
                            decoded_Bits = vitdec(rxSoft, trellis, traceback, ...
                                'trunc', 'unquant', puncturePattern);

                            decoded_Bits = decoded_Bits(1:payLoad_UncodedBits);

                        elseif strcmp(ChannelCoding, "Polar")

                            E = N_syms_perfram * M_bits;
                            nMax = 10;
                            listSize = 8;
                            iIL = false;
                            iBIL = true;

                            rxSoft = reshape(qamdemod(estdata, M_mod, 'gray', ...
                                'OutputType', 'llr', 'NoiseVariance', noiseVar), ...
                                N_bits_perfram, 1);

                            recPolar = nrRateRecoverPolar(rxSoft, K_crc, E, iBIL);
                            decoded_crc = nrPolarDecode(recPolar, K_crc, E, listSize, nMax, iIL, crcLen);
                            decoded_Bits = decoded_crc(1:payLoad_UncodedBits);

                        elseif strcmp(ChannelCoding, "None")

                            decoded_Bits = estBits;

                        else
                            error('Unknown ChannelCoding during decode.');
                        end

                    else
                        % ---------------- SP MODE DECODING ----------------
                        if strcmp(ChannelCoding, "None")

                            decoded_Bits = estBits;

                        elseif strcmp(ChannelCoding, "Polar")

                            E = N_syms_perfram * M_bits;
                            nMax = 10;
                            listSize = 8;
                            iIL = false;
                            iBIL = true;

                            rxSoft = reshape(qamdemod(softSymUnit, M_mod, 'gray', ...
                                'OutputType', 'llr', ...
                                'NoiseVariance', noiseVar / max(Pd, eps), ...
                                'UnitAveragePower', true), ...
                                N_bits_perfram, 1);

                            recPolar = nrRateRecoverPolar(rxSoft, K_crc, E, iBIL);
                            decoded_crc = nrPolarDecode(recPolar, K_crc, E, listSize, nMax, iIL, crcLen);
                            decoded_Bits = decoded_crc(1:payLoad_UncodedBits);

                        elseif strcmp(ChannelCoding, "Conv")

                            trellis = poly2trellis(7, [171 133]);

                            rxSoft = reshape(qamdemod(softSymUnit, M_mod, 'gray', ...
                                'OutputType', 'llr', ...
                                'NoiseVariance', noiseVar / max(Pd, eps), ...
                                'UnitAveragePower', true), ...
                                N_bits_perfram, 1);

                            traceback = 35;
                            decoded_Bits = vitdec(rxSoft, trellis, traceback, ...
                                'trunc', 'unquant', puncturePattern);

                            decoded_Bits = decoded_Bits(1:payLoad_UncodedBits);

                        else
                            error('SP mode currently supports only "None", "Polar", and "Conv".');
                        end
                    end

                    %% ============================================
                    % ERROR COUNT
                    % ============================================
                    errors = sum(xor(decoded_Bits, trans_info_bit));

                    err_ber(iesn0) = err_ber(iesn0) + errors;
                    avg_ber(iesn0) = err_ber(iesn0) / length(trans_info_bit) / frameCounter;

                    if any(decoded_Bits ~= trans_info_bit)
                        err_fer(iesn0) = err_fer(iesn0) + 1;
                    end
                    avg_fer(iesn0) = err_fer(iesn0) / frameCounter;

                    %% ============================================
                    % Display
                    % ============================================
                    clc
                    disp('####################################################################')
                    fprintf('RCP-OTFS-(N,M,QAM size) = [%d, %d, %d]\n', N, M, M_mod);
                    fprintf('ChannelEstimation = %s | Detection = %s | Coding = %s\n', ...
                        char(ChannelEstimation), char(Detection), char(ChannelCoding));
                    if strcmp(ChannelEstimation, "SP")
                        fprintf('SP params: Pd = %.2f, Pp = %.2f, numOuterIter = %d\n', Pd, Pp, numOuterIter);
                    end
                    fprintf('Current SNR = %d dB\n', SNR_dB(iesn0));
                    fprintf('Frames processed at this SNR = %d\n', frameCounter);
                    fprintf('Average BER = %.6e\n', avg_ber(iesn0));
                    fprintf('Average FER = %.6e\n', avg_fer(iesn0));
                    fprintf('Average NMSE = %.6e\n', avg_nmse(iesn0));

                end
            end
        end

        FER_curves{curveCount}  = avg_fer;
        BER_curves{curveCount}  = avg_ber;
        NMSE_curves{curveCount} = avg_nmse;
        curveLabels{curveCount} = sprintf('%s | %s | %s', ...
            char(ChannelCoding), char(Detection), char(ChannelEstimation));
    end
end

%% =========================
% Plots
% =========================
figure;
hold on;
for i = 1:curveCount
    cidx = mod(i-1, numel(colors)) + 1;
    midx = mod(i-1, numel(marker)) + 1;
    semilogy(SNR_dB, FER_curves{i}, 'LineWidth', 1.5, ...
        'Color', colors{cidx}, 'Marker', marker{midx});
end
grid on;
xlabel('SNR(dB)');
ylabel('FER');
title('RCP-OTFS: FER');
legend(curveLabels, 'Interpreter', 'none', 'Location', 'best');

figure;
hold on;
for i = 1:curveCount
    cidx = mod(i-1, numel(colors)) + 1;
    midx = mod(i-1, numel(marker)) + 1;
    semilogy(SNR_dB, BER_curves{i}, 'LineWidth', 1.5, ...
        'Color', colors{cidx}, 'Marker', marker{midx});
end
grid on;
xlabel('SNR(dB)');
ylabel('BER');
title('RCP-OTFS: BER');
legend(curveLabels, 'Interpreter', 'none', 'Location', 'best');

figure;
hold on;
for i = 1:curveCount
    cidx = mod(i-1, numel(colors)) + 1;
    midx = mod(i-1, numel(marker)) + 1;
    semilogy(SNR_dB, NMSE_curves{i}, 'LineWidth', 1.5, ...
        'Color', colors{cidx}, 'Marker', marker{midx});
end
grid on;
xlabel('SNR(dB)');
ylabel('Channel Operator NMSE');
title('RCP-OTFS: Channel NMSE');
legend(curveLabels, 'Interpreter', 'none', 'Location', 'best');

%% =========================
% Combined summary table
% =========================
SummaryTable = table(SNR_dB(:), 'VariableNames', {'SNR_dB'});

for i = 1:curveCount
    tag = matlab.lang.makeValidName(curveLabels{i});

    SummaryTable.(sprintf('BER_%s', tag))  = BER_curves{i}(:);
    SummaryTable.(sprintf('FER_%s', tag))  = FER_curves{i}(:);
    SummaryTable.(sprintf('NMSE_%s', tag)) = NMSE_curves{i}(:);
end

fprintf('\n============================================================\n');
fprintf('Combined Summary Table\n');
fprintf('============================================================\n');
disp(SummaryTable);

%% ============================================================
% Local helper functions for SP mode
% ============================================================

function [A, Gbasis] = build_sp_regression_matrix_rcp(N, M, l_i, k_i, s_ref)
    P = length(l_i);
    MN = N*M;

    A = zeros(MN, P);
    Gbasis = cell(P, 1);

    for p = 1:P
        Gbasis{p} = build_single_tap_Grcp(N, M, l_i(p), k_i(p));
        A(:, p) = Gbasis{p} * s_ref;
    end
end

function Gp = build_single_tap_Grcp(N, M, ell, k)
    MN = N*M;
    z = exp(1i*2*pi/MN);

    Gp = zeros(MN, MN);

    for q = 0:MN-1
        Gp(q+1, mod(q-ell, MN)+1) = z^(k * (q - ell));
    end
end

function [estBits, estdata, s_data_hat, softSymUnit] = sp_lmmse_detector_rcp( ...
    N, M, M_mod, noiseVar, r, G_hat, s_pilot, Pd, Modulation, lambda_x)

    MN = N*M;

    Fn = dftmtx(N);
    Fn = Fn ./ norm(Fn);

    Fm = dftmtx(M);
    Fm = Fm ./ norm(Fm);

    yData = r - G_hat * s_pilot;

    alpha = noiseVar / max(Pd, eps) + lambda_x;
    s_data_lin = (G_hat' * G_hat + alpha * eye(MN)) \ (G_hat' * yData);

    X_tilda_hat = reshape(s_data_lin, M, N);

    if(Modulation == 0)
        X_hat = X_tilda_hat * Fn;
    elseif(Modulation == 1)
        X_hat = Fm * X_tilda_hat;
    else
        error('Unknown Modulation option.');
    end

    softSymUnit = X_hat(:) / sqrt(Pd);

    idxHat = qamdemod(softSymUnit, M_mod, 'gray', ...
        'OutputType', 'integer', 'UnitAveragePower', true);

    estdata_unit = qammod(idxHat, M_mod, 'gray', ...
        'InputType', 'integer', 'UnitAveragePower', true);

    estdata = sqrt(Pd) * estdata_unit;

    estBits = reshape(qamdemod(estdata/sqrt(Pd), M_mod, 'gray', ...
        'OutputType', 'bit', 'UnitAveragePower', true), [], 1);

    Xhard = reshape(estdata, M, N);

    if(Modulation == 0)
        Xhard_tilda = Xhard * Fn';
    elseif(Modulation == 1)
        Xhard_tilda = Fm' * Xhard;
    else
        error('Unknown Modulation option.');
    end

    s_data_hat = reshape(Xhard_tilda, MN, 1);
end