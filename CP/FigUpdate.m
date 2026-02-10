%% Open existing .fig
h = openfig('LMMSEvsMRC_Perfect_SOFT.fig');   % Open saved figure
fig_num = h.Number;                           % Get figure number

% Get handles of objects
ax = findobj(h,'Type','Axes');                % Find axes handle(s)
ax = ax(1);                                   % Use first axes
lgd = findobj(h,'Type','Legend');             % Find legend handle (if any)
t   = get(ax,'Title');                        % Title handle
xlab= get(ax,'XLabel');                       % X-label handle
ylab= get(ax,'YLabel');                       % Y-label handle

% Apply formatting
if ~isempty(lgd)
    legend_proc(lgd);
end
plot_handles = findobj(ax,'Type','Line');     % All line objects in axes

plot_proc(h, plot_handles, '/home/sameer/OTFS/OTFS/SIMULATIONS/CP_OTFS', ...
           'LMMSEvsMRC_Perfect_SOFT', xlab, ylab, t);

%% Functions
function status = plot_proc(fig_handle, plot_handle, file_path, file_name, x_label, y_label, title_font)
    set(plot_handle, 'LineWidth', 2);
    set([x_label, y_label, title_font], 'FontSize', 16, 'Interpreter', 'latex');
    full_path = fullfile(file_path, file_name);
    ps_file = sprintf('%s.ps', full_path);
    pdf_file = sprintf('%s.pdf', full_path);

    print(fig_handle, '-dpsc2', ps_file);
    system(sprintf('ps2pdf %s.ps %s.pdf', full_path, full_path));
    system(sprintf('pdfcrop --margins "5 5 5 5" %s %s', pdf_file, pdf_file));

    status = 0;
end

function status = legend_proc(lgd)
    lgd.FontName = 'Arial';
    lgd.FontSize = 10;
    lgd.Interpreter = 'latex';
    lgd.Location = 'northeast';
    status = 0;
end
