function [est_bits,x_data] = LMMSE_detector(N,M,M_mod,noise_var,data_grid,r,G_rcp,Variant)
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix

Fm=dftmtx(M);  % Generate the DFT matrix
Fm=Fm./norm(Fm);
%% Initial assignments
%Number of symbols per frame
N_syms_perfram=sum(sum((data_grid>0)));
%Arranging the delay-Doppler grid symbols into an array
data_array=reshape(data_grid,1,N*M);
%finding position of data symbols in the array
[~,data_index]=find(data_array>0);
%number of bits per QAM symbol
M_bits=log2(M_mod);
%number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;
%received time domain blocks 

%%  LMMSE detection
   
    R=G_rcp'*G_rcp;
    sn_est=(R+noise_var.*eye(M*N))^(-1)*(G_rcp'*r);

X_tilda_est=reshape(sn_est,M,N);
%% detector output

switch Variant 
    case "OTFS"
        X_est=X_tilda_est*Fn;
    case "OFDM"
        X_est=Fm*X_tilda_est;
end

x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);
end
