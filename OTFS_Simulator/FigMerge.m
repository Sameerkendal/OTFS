% Open both figures invisibly
fig1 = openfig('/home/sameer/OTFS/OTFS/OTFS/OTFS_Simulator/Results_Figs/FER_OTFS_SPestimate_DiffCoding_0pt5.fig','invisible');
%fig2 =openfig('payload32_DiffCoding_perfectLMMSE_Rate0pt5_ofdm.fig','invisible');
fig2 = openfig('/home/sameer/OTFS/OTFS/OTFS/OTFS_Simulator/Results_Figs/FER_OTFS_Perfestimate_DiffCoding_0pt5.fig','invisible');

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
title(ax1,'Effect of Imperfect CSI over coded OTFS','FontSize',14,'FontWeight','bold')

%% -------- Edit Legend ----------
% Get all line handles in merged figure
lines = findall(ax1,'Type','line');



grid on
box on
