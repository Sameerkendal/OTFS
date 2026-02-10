function [H_t_f]=Generate_time_frequency_channel_CP(N,M,gs,L_set,length_cp)
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