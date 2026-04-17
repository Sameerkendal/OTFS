function [est_bits, x_data] = LMMSE_low_complexity(N,M,M_mod,noise_var,data_grid,r,G_rcp,mod)
%% Normalized DFT matrices
Fn = dftmtx(N);  Fn = Fn./norm(Fn);
Fm = dftmtx(M);  Fm = Fm./norm(Fm);

%% initilization
N_syms_perfram = sum(sum((data_grid>0)));
data_array     = reshape(data_grid,1,N*M);
[~,data_index] = find(data_array>0);

M_bits          = log2(M_mod);
N_bits_perfram  = N_syms_perfram*M_bits;

MN = M*N;

% Construct Psi = G*G' + sigma^2 I
Psi = G_rcp*G_rcp' + noise_var*eye(MN);

% Estimate delay spread (alpha) from sparsity of G_rcp
bandwidth = max(sum(abs(G_rcp)>0,2)); 
theta_blk = max(bandwidth-1,1);  

Q_sz = MN - theta_blk;

% Block partition
T_blk = sparse(Psi(1:Q_sz, 1:Q_sz));
B_blk = Psi(1:Q_sz, Q_sz+1:end);
S_blk = Psi(Q_sz+1:end, 1:Q_sz);
C_blk = Psi(Q_sz+1:end, Q_sz+1:end);

% LU decomposition (T block)
[L_T, U_T, P_T] = lu(T_blk);

% Compute E and V 
E_blk = L_T \ (P_T * B_blk);
V_blk = (U_T' \ S_blk')';

% Schur complement
FG_blk = C_blk - V_blk * E_blk;
[L_FG, U_FG, P_FG] = lu(FG_blk);

% Split received vector
r_top = r(1:Q_sz);
r_bot = r(Q_sz+1:end);

% Forward substitution 
y1 = L_T \ (P_T * r_top);
y2 = L_FG \ (P_FG * (r_bot - V_blk * y1));

% Backward substitution 
x2 = U_FG \ y2;
x1 = U_T \ (y1 - E_blk * x2);

% Final estimate 
sn_est = [x1; x2];

X_tilda_est = reshape(sn_est,M,N);

if(mod == 0)
    X_est = X_tilda_est*Fn;
else
    X_est = Fm*X_tilda_est;
end

x_est  = reshape(X_est,1,N*M);
x_data = x_est(data_index);

est_bits = reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'), ...
                   N_bits_perfram,1);

end