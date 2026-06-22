function [estBits, estdata] = MP_detector_RCP(r, N, M, M_mod, chan_coef, delay_taps, Doppler_taps,noiseVar)

MN = N*M;
Q  = M_mod;
P  = length(delay_taps);

Ni = 5;                                % No. of neighbors
nIter = 20;                            % No. of iterations

Fn = dftmtx(N);
Fn = Fn ./ norm(Fn);

Y_tilda = reshape(r, M, N);
Y_dd    = Y_tilda * Fn;
y       = reshape(Y_dd, MN, 1);

n_int = round(Doppler_taps);
kappa = Doppler_taps - n_int;          % residue in (-1/2, 1/2]

% Precompute Dirichlet table  D_N(q - kappa_i)  for q = -Ni..Ni
% Note: actual argument is (n_v - n_pp + k_i) = q + kappa_i  where q = (n_v - n_center)
% and n_center = n_pp - n_int_i. So beta = q + kappa_i.
alpha_tbl = zeros(2*Ni+1, P);
for i = 1:P
    for q = -Ni:Ni
        beta = q + kappa(i);
        if abs(beta - round(beta)) < 1e-12 && mod(round(beta), N) == 0
            % integer multiple of N -- only beta=0
            if abs(beta) < 1e-9
                alpha_tbl(q+Ni+1, i) = 1;     % D_N(0)/N = 1
            else
                alpha_tbl(q+Ni+1, i) = 0;
            end
        else
            alpha_tbl(q+Ni+1, i) = (exp(1j*2*pi*beta) - 1) / ...
                                   (N * (exp(1j*2*pi*beta/N) - 1));
        end
    end
end

% Build edge list 
Emax = MN * P * (2*Ni + 1);
edges_d = zeros(1, Emax);
edges_c = zeros(1, Emax);
edges_h = zeros(1, Emax);
E = 0;

for mp = 0:M-1
    for npp = 0:N-1
        d = mp + npp*M + 1;

        for i = 1:P
            l_i = delay_taps(i);
            k_i = Doppler_taps(i);
            m_v = mod(mp - l_i, M);
            s_borrow = double(mp < l_i);

            delay_phase = exp(1j*2*pi * k_i * (mp - l_i) / MN);

            n_center = mod(npp - n_int(i), N);

            for q = -Ni:Ni
                alpha = alpha_tbl(q+Ni+1, i);
                if abs(alpha) < 1e-12, continue; end

                n_v = mod(n_center + q, N);
                c = m_v + n_v*M + 1;

                borrow_phase = exp(-1j*2*pi*s_borrow*n_v/N);

                E = E + 1;
                edges_d(E) = d;
                edges_c(E) = c;
                edges_h(E) = chan_coef(i) * delay_phase * borrow_phase * alpha;
            end
        end
    end
end

edges_d = edges_d(1:E);
edges_c = edges_c(1:E);
edges_h = edges_h(1:E);

% Adjacency
edges_at_d = cell(MN, 1);
edges_at_c = cell(MN, 1);
for e = 1:E
    edges_at_d{edges_d(e)}(end+1) = e;
    edges_at_c{edges_c(e)}(end+1) = e;
end

A = qammod((0:Q-1).', Q, 'gray');

% MP
msg = ones(E, Q) / Q;

Delta     = 0.6;
gamma     = 0.01;
eta_best  = -1;
best_xhat = zeros(MN, 1);

for iter = 1:nIter

    mu_total  = zeros(MN, 1);
    var_total = zeros(MN, 1);

    for d = 1:MN
        elist = edges_at_d{d};
        for kk = 1:length(elist)
            e = elist(kk);
            h = edges_h(e);
            p = msg(e, :);
            Ep  = p * A;
            Ep2 = p * abs(A).^2;
            mu_total(d)  = mu_total(d)  + h * Ep;
            var_total(d) = var_total(d) + abs(h)^2 * (Ep2 - abs(Ep)^2);
        end
    end

    log_marg_all = zeros(MN, Q);

    for c = 1:MN
        elist = edges_at_c{c};
        Sc = length(elist);
        log_xi = zeros(Sc, Q);

        for jj = 1:Sc
            e = elist(jj);
            d = edges_d(e);
            h = edges_h(e);
            p = msg(e, :);
            Ep  = p * A;
            Ep2 = p * abs(A).^2;

            mu_dc  = mu_total(d)  - h * Ep;
            var_dc = var_total(d) - abs(h)^2 * (Ep2 - abs(Ep)^2) + noiseVar;
            var_dc = max(real(var_dc), 1e-10);

            resid = y(d) - mu_dc - h * A.';
            log_xi(jj, :) = -abs(resid).^2 / var_dc;
        end

        log_sum = sum(log_xi, 1);
        log_marg_all(c, :) = log_sum;

        for jj = 1:Sc
            e = elist(jj);
            log_p = log_sum - log_xi(jj, :);
            log_p = log_p - max(log_p);
            p_new = exp(log_p);
            p_new = p_new / (sum(p_new) + eps);
            msg(e, :) = Delta * p_new + (1-Delta) * msg(e, :);
        end
    end

    log_marg_all = log_marg_all - max(log_marg_all, [], 2);
    marg = exp(log_marg_all);
    marg = marg ./ (sum(marg, 2) + eps);

    [maxp, idx] = max(marg, [], 2);
    xhat_iter = idx - 1;
    eta = mean(maxp >= 1 - gamma);

    if eta > eta_best
        eta_best  = eta;
        best_xhat = xhat_iter;
    end
    if eta >= 1, break; end
end

estdata = qammod(best_xhat, Q, 'gray');
estBits = reshape( ...
    qamdemod(estdata, Q, 'gray', 'OutputType','bit'), ...
    [], 1);

end