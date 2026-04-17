function [est_bits,x_data] = Block_LMMSE_detector(N,M,M_mod,noise_var,data_grid,r,gs,L_set,length_cp,U,Modulation)
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
Fm=dftmtx(M);  % Generate the DFT matrix
Fm=Fm./norm(Fm);  % normalize the DFT matrix
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
sn_block_est=zeros(M,N);
Gn=zeros(M*U,M*U);
%% block-wise LMMSE detection
for n=1:N    
    for m=1:M*U
        for l=L_set
                Gn(m,mod(m-l-1,M*U)+1)=gs(l+1,(m-1)+ (n-1)*(M+length_cp)*U+1);           
        end
    end
    rn=r((n-1)*M*U+1:n*M*U);    
    Rn=Gn'*Gn;
    sn_block_est(:,n)=resample((Rn+noise_var.*eye(M*U))^(-1)*(Gn'*rn),1,U);
end
X_tilda_est=sn_block_est;

%% detector output
if(strcmp(Modulation,"OTFS"))
    X_est=X_tilda_est*Fn;
elseif(strcmp(Modulation,"OFDM"))
    X_est=Fm*X_tilda_est;
end

x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

end
