function [H_t_f]=Generate_time_frequency_channel(N,M,G)

    Fmn=dftmtx(M*N);
    Fmn=Fmn./norm(Fmn);
    
    H_t_f = Fmn*G*Fmn';
end