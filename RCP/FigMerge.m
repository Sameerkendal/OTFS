% Open both figures invisibly
fig1 = openfig('/home/sameer/OTFS/OTFS/SIMULATIONS/RCP_OTFS/DiffCoding_otfs.fig','invisible');
%fig2 =openfig('payload32_DiffCoding_perfectLMMSE_Rate0pt5_ofdm.fig','invisible');
fig2 = openfig('/home/sameer/OTFS/OTFS/SIMULATIONS/CP_OTFS/Diffcoding_ofdm.fig','invisible');

% Get axes
ax1 = findall(fig1,'Type','axes');
ax2 = findall(fig2,'Type','axes');

% Merge plots from fig2 into fig1
h2 = copyobj(allchild(ax2), ax1);

% Bring figure to front
figure(fig1)
hold(ax1,'on')

%% -------- Edit Line Styles (optional) ----------
% Example: make fig2 lines dotted
set(h2,'LineStyle',':','LineWidth',1.5)

%% -------- Edit Title ----------
title(ax1,'Effect of Error Correction Codes','FontSize',14,'FontWeight','bold')

%% -------- Edit Legend ----------
% Get all line handles in merged figure
lines = findall(ax1,'Type','line');

% Define legend entries (order matters!)
% legend(ax1, lines, ...
%     {'OFDM : 5Tap','OFDM : 3Tap','OFDM : 2Tap','OTFS: 5Tap','OTFS : 3Tap','OTFS : 2Tap','OFDM : Convolutional','OTFS : No coding','OTFS : Polar','OTFS : Convolutional'}, ...
%     'Location','southwest')

legend(ax1, lines, ...
    {'OFDM : None','OFDM : Polar','OFDM : Convolutional','OTFS : None','OTFS : Polar','OTFS : Convolutional','OTFS: 5Tap','OTFS : 3Tap','OTFS : 2Tap','OFDM : Convolutional','OTFS : No coding','OTFS : Polar','OTFS : Convolutional'}, ...
    'Location','southwest')

grid on
box on
