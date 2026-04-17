function [H_t_f]=Generate_time_frequency_channel(N,M,G,gs,L_set,length_cp,Variant)

    switch Variant
        case 'RCP'
            Fmn=dftmtx(M*N);
            Fmn=Fmn./norm(Fmn);
            
            H_t_f = Fmn*G*Fmn';

        case 'CP'
            H_t_f=zeros(N,M); % Time-frequency single tap channel matrix
            Fm=dftmtx(M);
            Fm=Fm./norm(Fm);
            
            for n=1:N
                Gn=zeros(M,M);
                for m=1:M
                    for l=L_set
                        Gn(m,mod(m-l-1,M)+1)=gs(l+1,(m-1)+ (n-1)*(M+length_cp)+1);
                    end
                end
                H_t_f(n,1:M)=diag(Fm*Gn*Fm').';
            end

    end
end