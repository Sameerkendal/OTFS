function [est_bits,ite,x_data] = MRC_delay_time_detector(N,M,M_data,M_mod,no,data_grid,r,H_tf,nu_ml_tilda,L_set,omega,decision,init_estimate,n_ite_MRC)

    Fn=dftmtx(N);  
    Fn=Fn./norm(Fn);  
    
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
    L_set=L_set+1; 

%% initial estimate using single tap TF equalizer
    if(init_estimate==1)
        Y_tf=fft(Y_tilda).';                                                                                                       
        X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+no);                                                             
        X_est = ifft(X_tf.')*Fn;                                                                                               
        X_est=qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray');
        X_est=X_est.*data_grid;
        X_tilda_est=X_est*Fn';
    else
        X_est=zeros(M,N);
        X_tilda_est=X_est*Fn';
    end
    x_m=X_est.';
    x_m_tilda=X_tilda_est.';


%% MRC detector    %% Algorithm 2 in [R1] (or Algotithm 3 in Chapter 6, [R2])

    d_m_tilda=zeros(N,M);
    y_m_tilda=reshape(r,M,N).';
    delta_y_m_tilda=y_m_tilda;
    for m=1:M_prime   
        for l=L_set
            d_m_tilda(:,m)=d_m_tilda(:,m)+abs(nu_ml_tilda(:,(mod(m+(l-1)-1,M)+1),l).^2);                                                             
        end
    end
    for m = 1:M
        for l = L_set
            idx_x = mod(m - (l-1) - 1, M) + 1;  
            delta_y_m_tilda(:, m) = delta_y_m_tilda(:, m) - nu_ml_tilda(:, m, l) .* x_m_tilda(:, idx_x);
        end
    end
    x_m_tilda_old=x_m_tilda;
    c_m_tilda=x_m_tilda;

%% iterative computation
    for ite=1:(n_ite_MRC)                                                                                                                 
        delta_g_m_tilda=zeros(N,M);
        for m=1:M_prime         
            for l=L_set
                 
                idx_delta_y = mod(m + (l-1) - 1, M) + 1;
                idx_nu = mod(m + (l-1) - 1, M) + 1; 
                delta_g_m_tilda(:, m) = delta_g_m_tilda(:, m) + conj(nu_ml_tilda(:, idx_nu, l)) .* delta_y_m_tilda(:, idx_delta_y);
            end
            c_m_tilda(:,m)=x_m_tilda_old(:,m)+delta_g_m_tilda(:,m)./d_m_tilda(:,m);                                                     
            if(decision==1)
                x_m(:,m)=qammod(qamdemod(Fn*(c_m_tilda(:,m)),M_mod,'gray'),M_mod,'gray');                                               
                x_m_tilda(:,m)=(1-omega)*c_m_tilda(:,m)+omega*Fn'*x_m(:,m);
            else
                x_m_tilda(:,m)=c_m_tilda(:,m);
            end
            for l=L_set                                                                                                                
                idx_target = mod(m + (l-1) - 1, M) + 1;  
                delta_y_m_tilda(:, idx_target) = delta_y_m_tilda(:, idx_target) - nu_ml_tilda(:, idx_target, l) .* (x_m_tilda(:, m) - x_m_tilda_old(:, m));
            end                                                                                                                        
            x_m_tilda_old(:,m)=x_m_tilda(:,m);
        end
           
        %% convergence criteria
        error(ite)=norm(delta_y_m_tilda);
        if(ite>1)
            if(error(ite)>=error(ite-1))                                                                                                 
                break;
            end
        end   
    end
        %%
    if(n_ite_MRC==0)
        ite=0;
    end
    %% detector output bits
    X_est=(Fn*x_m_tilda).';
    x_est=reshape(X_est,1,N*M);
    x_data=x_est(data_index);
    est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);                                                % Line 17 of Algorithm 2 in [R1]

end
