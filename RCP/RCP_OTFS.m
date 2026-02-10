close all
clear all

%%
colors = ["r", "g",  "c", "m", "y", "k", 'r', 'g',"b", 'b', ...
          'c', 'm', 'y', 'k', 'r', 'g', 'b', 'c', 'm', 'y'];

marker =  {"o","+","*","square","diamond","^","x","pentagram","v",">","."};

%% OTFS parameters%%%%%%%%%%

N = 64;
M = 4;
M_mod = 4;
M_bits = log2(M_mod);
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
BandWidth = 2*10^6;
TdMax = 0.5*10^-6;
data_grid=ones(M,N);


N_syms_perfram = sum(sum(data_grid));
N_bits_perfram = N_syms_perfram*M_bits;


car_fre = 18*10^9;
delta_f= BandWidth/M; 
T = 1/delta_f; 
payLoad_UncodedBits = 384;
codeRate = 0.75;
SNR_dB = 0:5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;


N_fram = 10;
seedGen =1:50:5000;
d=0;
%% Normalized DFT matrix
Fn=dftmtx(N);  
Fn=Fn./norm(Fn);  
Fm=dftmtx(M); 
Fm=Fm./norm(Fm); 

ChannelEstimation =["Perfect"];
Modulation = 0;
%"OTFS";  
%ChannelCoding = "Polar";
%ChannelCoding = ["None"];
for ChannelCoding=["Conv","Polar"]
    for Detection = ["LMMSE"]
        d=d+1;
        
        current_frame_number=zeros(1,length(SNR_dB));
        estBits=zeros(N_bits_perfram,1);
        err_ber = zeros(1,length(SNR_dB));
        avg_ber = zeros(1,length(SNR_dB));
        err_fer = zeros(1,length(SNR_dB));
        avg_fer = zeros(1,length(SNR_dB));

    
        for iesn0 = 1:length(SNR_dB)
            for a =1:length(seedGen)
                seed = seedGen(a);
                for ifram = 1:N_fram
                    current_frame_number(iesn0)= ifram;

                    %% FRAME GENERATION %%%%%
                    trans_info_bit = randi([0,1],payLoad_UncodedBits,1);  
                    
                    %% ENCODING %%%%%%
                    outLen = 512;
                    
                    

                    if(strcmp(ChannelCoding,"Conv"))
                         padLen = (outLen * codeRate) - payLoad_UncodedBits;

                            if(padLen ~= 0)
                                 padBits = zeros(padLen, 1);
                                 mssgBits = [trans_info_bit; padBits];
                            else
                                mssgBits = trans_info_bit;
                            
                            end      
                        if(codeRate == 0.5)
                        %%% RATE = 1/2 %%%%
                            puncturePattern = [1,1 ];
                        elseif(codeRate == 0.75)
                        %%% RATE = 3/4 %%%%
                            puncturePattern = [1 1 1 0 0 1];   
                        end

                            trellis = poly2trellis(7, [171 133]);
                            
                               
                            encodedBits = convenc(mssgBits, trellis,puncturePattern);
                            if length(encodedBits) >= outLen
                                codedBits = encodedBits(1:outLen);       
                            else
                                reps = ceil(outLen / length(encodedBits));
                                codedBits = repmat(encodedBits, reps, 1);
                                codedBits = codedBits(1:outLen);         
                            end
                    
                            mappedBits = reshape(codedBits, M_bits, []);     
                            data = qammod(mappedBits, M_mod, 'gray', 'InputType', 'bit');
                    
                           
                            X = Generate_2D_data_grid(N, M, data, data_grid);
                            
                    elseif(strcmp(ChannelCoding,"Polar"))
                            mssgBits = trans_info_bit;
                            % padLen = (outLen * codeRate) - payLoad_UncodedBits;
                            % 
                            % if(padLen ~= 0)
                            %      padBits = zeros(padLen, 1);
                            %      mssgBits = [trans_info_bit; padBits];
                            % else
                            %     mssgBits = trans_info_bit;
                            % 
                            % end
                            K = length(mssgBits);                    
                            crcLen = 11;  
                            msgCRC = nrCRCEncode(mssgBits, '11');   
                            K_crc = length(msgCRC);
                            E = N_syms_perfram * M_bits;      
                            nMax = 9;                         
                            iIL = false;                       
                        
                            enc = nrPolarEncode(msgCRC, E, nMax, iIL);  
                        
                            modIn = nrRateMatchPolar(enc, length(msgCRC), E, true);  
                          
                            mappedBits = reshape(modIn, M_bits, []);
                            data = qammod(mappedBits, M_mod, 'gray', 'InputType', 'bit');
                        
                            
                            X = Generate_2D_data_grid(N, M, data, data_grid);
                            
                    elseif(strcmp(ChannelCoding,"None"))
                            trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
                            data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
                            X = Generate_2D_data_grid(N,M,data,data_grid);
                    end
                    
                    %% Modulation%%%%
                    if(Modulation == 0)   %%OTFS
                        X_tilda=X*Fn'; 
                    elseif(Modulation==1)  %%OFDM
                        X_tilda = Fm'*X;
                    end
                    s = reshape(X_tilda,N*M,1);   
                    
                    %% CHANNEL CONSTRUCTION
    
                    max_speed=12000;  
                    [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed,seed,1);
                      
                    L_set=unique(delay_taps);       
                    
                    Lmax = max(delay_taps);
                    
                    gs=zeros(Lmax+1,N*(M));
                    z=exp(1i*2*pi/(N*(M)));
                    l_i = delay_taps;
                    g_i=chan_coef;
                    k_i=Doppler_taps;
                    
                    for q=0:N*(M)-1
                        for i=1:taps
                            gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z^(k_i(i)*(q-l_i(i)));
                        end
                    end
    
                    G_rcp=zeros(N*M,N*M);
                    for q=0:N*M-1
                        for ell=0:Lmax
                               G_rcp(q+1,mod(q-ell,N*M)+1)=gs(ell+1,q+1);
                        end
                    end
                   
                    
                     r=G_rcp*s;
                     
                    %% channel output%%%%%   
    
                     noise= sqrt(sigma_2(iesn0)/2)*(randn(size(r)) + 1i*randn(size(r)));
                     r = r+noise;
                     
                     %% OTFS DEMODULATION %%%%%%
    
                     Y_tilda=reshape(r,M,N); 
                     if(Modulation==0)   %%OTFS
                        Y = Y_tilda*Fn; 
                    elseif(Modulation==1) %%OFDM
                        Y = Fm*Y_tilda;
                    end
                     
                    if(strcmp(ChannelEstimation,"Estimator")) 
%% PILOT TRANSMISSION


                            xp = 6;
                            pilot_grid = zeros(M,N);
                            pilot_grid(2,32) = 1;
                            pilot_positions = find(pilot_grid);
            
                            N_pilot_syms = sum(sum(pilot_grid));
                            pilot_data =xp* ones(N_pilot_syms, 1);
                            X_pilot = Generate_2D_data_grid(N,M,pilot_data,pilot_grid);
            
                            X_pilot_tilda=X_pilot*Fn';                     
                            s_pilot = reshape(X_pilot_tilda,N*M,1);
            
                            r_pilot=G_rcp*s_pilot;
                            noise_pilot = sqrt(sigma_2(iesn0)/2)*(randn(size(s_pilot)) + 1i*randn(size(s_pilot)));
                            r_pilot=r_pilot+noise_pilot;
            
                            Y_pilot_tilda=reshape(r_pilot,M,N);
                            Y_pilot = Y_pilot_tilda*Fn;

%% CHANNEL ESTIMATION %%%%%

                            max_delay_taps = max(delay_taps);
                            max_doppler_taps = 31;
        
                            [pilot_m, pilot_n] = ind2sub([M, N], pilot_positions(1));
        
        
                            k_est_range = max(1, pilot_n - max_doppler_taps) : min(N, pilot_n + max_doppler_taps);
                            l_est_range = max(1, pilot_m) : min(M, pilot_m + max_delay_taps);
        
                            threshold = 5*sqrt(sigma_2(iesn0));  
        
        
        
                            l_i = [];
                            k_i = [];
                            g_i = [];
                            count = 1;
        
                            for k = k_est_range
                                for l = l_est_range
                                    y_received = Y_pilot(l, k);
        
                                    if abs(y_received) > threshold
                                        k_rel = k - pilot_n;
                                        l_rel = l - pilot_m;
                                        phase_compensation = exp(-1j * 2 * pi*(1/(M*N)) * k_rel * (pilot_m-1));
                                        %phase_compensation=1;
                                        l_i(count) = l_rel;
                                        k_i(count) = k_rel;
                                        g_i(count) = (y_received / xp) * phase_compensation;  % This is the COMBINED gain at this (l,k)
                                        count = count + 1; 
                                    end
                                end
                            end 

%% Channel Reconstruction

                        taps = length(l_i);
                        z = exp(1i*2*pi/(N*(M)));
                        delay_spread = max_delay_taps;
                        gs_est = zeros(delay_spread+1, N*(M));
        
                        for q = 0:(N*(M)-1)
                            for i = 1:taps
                                gs_est(l_i(i)+1, q+1) = gs_est(l_i(i)+1, q+1) + g_i(i) * z^(k_i(i)*(q-l_i(i)));
                            end
                        end
                        
                        G_rcp_est = zeros(N*M,N*M);
                            for q = 0:N*M-1
                                for ell = 0:Lmax
                                    G_rcp_est(q+1,mod(q-ell,N*M)+1) = gs_est(ell+1,q+1);
                                end
                            end
        
                   end
         

                    if(strcmp(ChannelEstimation,"Noisy"))
                   %% NOISY CHANNEL MODEL %%%
    
                            EstmationGain_dB = 10;
                            pilotSNR_dB = iesn0 + EstmationGain_dB;
                            pilotSNR = 10^(pilotSNR_dB/10);
                            Ppilot = 1;
                            Npilots = 1;
                            sigma_e2 = (Ppilot / pilotSNR) / (abs(Ppilot)^2 * Npilots);
                            
                            e_i = sqrt(sigma_e2/2) * (randn(size(g_i)) + 1i*randn(size(g_i)));
                            g_i_hat = g_i + e_i;
                            
                            gs_hat = zeros(Lmax+1, N*(M));
                            z = exp(1i*2*pi/(N*(M)));
                            
                            for q = 0:N*(M)-1
                                for i = 1:taps
                                    gs_hat(l_i(i)+1, q+1) = gs_hat(l_i(i)+1, q+1) + g_i_hat(i) * z^(k_i(i)*(q - l_i(i)));
                                end
                            end

                            G_rcp_hat = zeros(N*M,N*M);
                            for q = 0:N*M-1
                                for ell = 0:Lmax
                                    G_rcp_hat(q+1,mod(q-ell,N*M)+1) = gs_hat(ell+1,q+1);
                                end
                            end
                    end

                     
%% Channel Vector Estimation %%%%%
    
                    if(strcmp(ChannelEstimation,"Perfect"))
                        gs = gs;
                        G_rcp = G_rcp;
                    elseif(strcmp(ChannelEstimation,"Estimator"))
                        gs = gs_est;
                        G_rcp = G_rcp_est;
                    elseif(strcmp(ChannelEstimation,"Noisy"))
                        gs = gs_hat;
                        G_rcp = G_rcp_hat;
                    end
                   [H,H_tilda,P]= Gen_DD_and_DT_channel_matrices(N,M,G_rcp,Fn);
                   [nu_ml_tilda,~,K_ml]=Gen_DT_and_DD_channel_vectors(N,M,Lmax,gs);   %DELAY TIME CHANNEL
                   [H_tf]=Generate_time_frequency_channel(N,M,G_rcp);
                     
                     
                   
                   %% DETECTION
                    if(strcmp(Detection,"MRC"))
                        n_ite_MRC=100; 
                        omega=1;
                        if(M_mod==64)
                            omega=0.25;     
                        end
                        decision=0;         
                        init_estimate=0;    
                        [estBits,det_iters_MRC,estdata] = MRC_delaydoppler(N,M,M,M_mod,sigma_2(iesn0),data_grid,r,H_tf,L_set,omega,decision,init_estimate,n_ite_MRC,K_ml,Y);
                    elseif(strcmp(Detection,"LMMSE"))
                        [estBits,estdata]= LMMSE_detector(N,M,M_mod,sigma_2(iesn0),data_grid,r,G_rcp,Modulation);
                       
                    end

                    %% DECODING
                    if(strcmp(ChannelCoding,"Conv"))
                    
                            trellis = poly2trellis(7, [171 133]);   
                            rxSoft = reshape(qamdemod(estdata, M_mod, 'gray','OutputType', 'llr','NoiseVariance', sigma_2(iesn0)),N_bits_perfram, 1);
                            
                            traceback = 35;   
                            
                                 decoded_Bits = vitdec(rxSoft, ...
                                                          trellis, ...
                                                          traceback, ...
                                                          'trunc', ...
                                                          'unquant', ...
                                                          puncturePattern);

                                    decoded_Bits = decoded_Bits(1:payLoad_UncodedBits);
                           
                           
                    
                    elseif(strcmp(ChannelCoding,"Polar"))

                        E       = N_syms_perfram * M_bits;
                        nMax    = 10;
                        listSize = 8;
                        iIL      = false;    
                        
                        rxSoft = reshape(qamdemod(estdata, M_mod, 'gray','OutputType', 'llr','NoiseVariance', sigma_2(iesn0)),N_bits_perfram, 1);
                        iBIL = true;
                        recPolar = nrRateRecoverPolar(rxSoft, K_crc, E, iBIL);
                        
                        decoded_crc = nrPolarDecode(recPolar, K_crc, E, listSize, nMax, iIL, crcLen);
                        
                        decoded_Bits = decoded_crc(1:payLoad_UncodedBits);

                    elseif(strcmp(ChannelCoding,"None"))
                    
                            decoded_Bits = estBits;
                    
                    end
                    
                    %% errors count %%
    
                    errors = sum(xor(decoded_Bits,trans_info_bit));                
                    err_ber(iesn0) = err_ber(iesn0) + errors;        
                    avg_ber(iesn0)=err_ber(iesn0).'/length(trans_info_bit)/(ifram+((a-1)*N_fram));
                    

                     %% Frame error check
                    if any(decoded_Bits ~= trans_info_bit)
                        err_fer(iesn0) = err_fer(iesn0) + 1;   
                    end
                    
                    
                    avg_fer(iesn0) = err_fer(iesn0) / (ifram + (a-1)*N_fram);

                    %% DISP error performance details   
    
                    clc
                    disp('####################################################################') 
                    fprintf('RCP-OTFS-(N,M,QAM size)');disp([N,M,M_mod]);
                    fprintf('RCP-OTFS - TWO TAP ')
                    display(SNR_dB,'SNR (dB)');
                    display(current_frame_number,'Number of frames');
                    display(avg_fer,'Average FER');                                
                end
            end
        end
        if(Modulation == 0)
            semilogy(SNR_dB,avg_fer,LineWidth=1.5,LineStyle='-',Color=colors{d},Marker=marker{d+1});
            hold on;
        else
            semilogy(SNR_dB,avg_fer,LineWidth=2,LineStyle=':',Color=colors{d},Marker=marker{d});
            hold on;
        end
    end
end
title("Performance with Different Coding schemes (Packet size : 48 bytes)")
%legend("OTFS : Mach 1","OTFS : Mach 2","OTFS : Mach 4","OTFS : Mach 8","OTFS : Mach 10")

%legend("CodeRate = 0.5 ","CodeRate = 0.75","None","LMMSE")
%legend("Payload = 240","Payload = 256","Payload = 288","Payload = 320","Payload = 352","Payload = 384")
legend("Convolutional","Polar","None")
grid on
xlabel('SNR(dB)')
ylabel('FER')


