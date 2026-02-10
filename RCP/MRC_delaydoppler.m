
function [est_bits,ite,x_data] = MRC_delaydoppler(N,M,M_data,M_mod,no,data_grid,r,H_tf,L_set,omega,decision,init_estimate,n_ite_MRC,K_ml,Y)
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
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
%received delay-time samples
Y_tilda=reshape(r,M,N);
M_prime=M_data;
L_set=L_set+1; % Since matlab cannot work with zero indices
%% initial estimate using single tap TF equalizer
if(init_estimate==1)
    Y_tf=fft(Y_tilda).'; % delay-time to frequency-time domain                                                                      % equation (63) in [R1]
    X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+no); % single tap equalizer                                                            % equation (64) in [R1]
    X_est = ifft(X_tf.')*Fn; % SFFT                                                                                                 % equation (65) in [R1]
    X_est=qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray');
    X_est=X_est.*data_grid;
else
    X_est=zeros(M,N);
end
x_m_hat=X_est.';

%% MRC detector    %% Algorithm 1 in [R1]
%% initial computation
D_m=zeros(N,N,M);
D_m_inv=zeros(N,N,M);
y_m=Y.';
% RNPI_error=y_m;
for m=1:M_prime
    for l=L_set
        D_m(:,:,m)=D_m(:,:,m)+K_ml(:,:,mod(m+(l-1)-1,M)+1,l)'*K_ml(:,:,mod(m+(l-1)-1,M)+1,l);                                                             % equation (47) in [R1]
    end
    D_m_inv(:,:,m)=inv(D_m(:,:,m));
end
c_m=x_m_hat;
b_m_l=zeros(N,M,length(L_set));

%% iterative computation
for ite=1:n_ite_MRC
    g_m=zeros(N,M);
    RNPI_error=y_m;
    for m=1:M_prime
        for l=L_set
            b_m_l(:,m,l)=y_m(:,mod(m+(l-1)-1,M)+1);
            for p=L_set
                if(l~=p)
                    
                        b_m_l(:,m,l)=b_m_l(:,m,l)-K_ml(:,:,mod(m+(l-1)-1,M)+1,p)*x_m_hat(:,mod(m+(l-p)-1,M)+1);                                % Line 5 of Algorithm 1 in [R1] 
                    
                end
            end
            g_m(:,m)=g_m(:,m)+K_ml(:,:,mod(m+(l-1)-1,M)+1,l)'*b_m_l(:,m,l);                                                       % Line 7 of Algorithm 1 in [R1] 
        end
        c_m(:,m)=D_m_inv(:,:,m)*g_m(:,m);                                                                              % Line 8 of Algorithm 1 in [R1] 
        if(decision==1)
            x_m_hat(:,m)=(1-omega)*c_m(:,m)+omega*qammod(qamdemod((c_m(:,m)),M_mod,'gray'),M_mod,'gray');              % Line 9 of Algorithm 1 in [R1] 
        else
            x_m_hat(:,m)=c_m(:,m);
        end
        
        for l=L_set
            
                RNPI_error(:,m)=RNPI_error(:,m)-reshape(K_ml(:,:,m,l),N,N)*x_m_hat(:,mod(m-(l-1)-1,M)+1);                         % equation (50) in [R1] - residual interference
            
        end
    end
    %% convergence criteria
    error(ite)=norm(RNPI_error);
    if(ite>1)
        if(error(ite)>=error(ite-1))
            break;
        end
    end
end
if(n_ite_MRC==0)
    ite=0;
end
%% detector output bits
X_est=x_m_hat.';
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

end
