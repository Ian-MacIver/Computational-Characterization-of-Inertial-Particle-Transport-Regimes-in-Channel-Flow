function plot_turbulent_particle_results_Final()

clearvars;
clc;

%% USER SETTINGS 
USE_ENSEMBLE    = true;
SINGLE_SNAPSHOT = 1;
ENSEMBLE_FILE   = 'turbulent_ensemble_results.mat';
RESULT_PREFIX   = 'turbulent_particle_final_results_';

%% Global plot defaults 
set(0, 'DefaultAxesFontSize',    13);
set(0, 'DefaultAxesFontName',    'Arial');
set(0, 'DefaultLineLineWidth',   1.6);
set(0, 'DefaultFigureColor',     'w');
set(0, 'DefaultAxesGridAlpha',   0.25);
set(0, 'DefaultAxesGridLineStyle', ':');

fig_size = [0, 0, 820, 580];

%% Load data 
if USE_ENSEMBLE
    if ~isfile(ENSEMBLE_FILE)
        error(['Could not find %s.\n' ...
               'Run compute_ensemble_statistics.m first.'], ENSEMBLE_FILE);
    end

    E = load(ENSEMBLE_FILE);

    tau_p_list          = E.tau_p_list;
    v0_list             = E.v0_list;
    St_list             = E.St_list;
    v0_norm             = E.v0_norm;
    params              = E.ref_params;

    ymax_slice          = E.ensemble_mean_ymax;
    ymax_std_slice      = E.ensemble_std_ymax;

    outcome_map         = E.outcome_stack;
    t_res_map           = E.t_res_stack;
    bounce_at_stick_map = E.bounce_stack;

    alpha_partial       = E.alpha_partial_ensemble;
    C_partial           = E.C_partial_ensemble;
    alpha_mid           = E.alpha_mid_ensemble;
    C_mid               = E.C_mid_ensemble;

    running_mean_alpha_mid     = E.running_mean_alpha_mid;
    running_mean_alpha_partial = E.running_mean_alpha_partial;
    converged_at               = E.converged_at;
    N_SNAPSHOTS                = E.N_SNAPSHOTS;

    % Trajectories may or may not exist in the ensemble file
    has_trajectories = isfield(E, 'trajectories') && ~isempty(E.trajectories);
    if has_trajectories
        trajectories = E.trajectories;
    end

    mode_label = 'Ensemble (11 snapshots)';
    fprintf('Loaded ensemble results from %s\n', ENSEMBLE_FILE);
    fprintf('Ensemble alpha_mid = %.4f,  alpha_partial = %.4f\n', ...
        alpha_mid, alpha_partial);
    fprintf('Convergence confirmed at snapshot %d of %d\n', ...
        converged_at, N_SNAPSHOTS);

else
    results_file = sprintf('%s%d.mat', RESULT_PREFIX, SINGLE_SNAPSHOT);
    if ~isfile(results_file)
        error('Could not find %s. Run the sweep for snapshot %d first.', ...
              results_file, SINGLE_SNAPSHOT);
    end

    S = load(results_file);

    tau_p_list          = S.tau_p_list;
    v0_list             = S.v0_list;
    St_list             = S.St_list;
    params              = S.params;
    outcome_map         = S.outcome_map;
    bounce_at_stick_map = S.bounce_at_stick_map;
    t_res_map           = S.t_res_map;

    Nx0 = numel(S.x0_list);
    ix  = ceil(Nx0 / 2);
    ymax_slice     = squeeze(S.y_max_map(:,:,ix));
    ymax_std_slice = [];

    v0_norm = v0_list ./ params.ref.Uref;

    [alpha_partial, C_partial] = fit_boundary_local( ...
        St_list, v0_norm, ymax_slice, -0.5, 4);
    [alpha_mid, C_mid] = fit_boundary_local( ...
        St_list, v0_norm, ymax_slice,  0.0, 4);

    running_mean_alpha_mid     = [];
    running_mean_alpha_partial = [];
    converged_at               = [];
    N_SNAPSHOTS                = 1;

    has_trajectories = isfield(S, 'trajectories') && ~isempty(S.trajectories);
    if has_trajectories
        trajectories = S.trajectories;
    end

    mode_label = sprintf('Single snapshot (t = %d)', SINGLE_SNAPSHOT);
    fprintf('Loaded single snapshot %d from %s\n', SINGLE_SNAPSHOT, results_file);
end

geom = params.geom;
ref  = params.ref;

Ntau = numel(tau_p_list);
Nv0  = numel(v0_list);

fprintf('\nMode: %s\n', mode_label);
fprintf('Parameter space: %d x %d  (%d total cases)\n', Ntau, Nv0, Ntau*Nv0);

%% Output folder
if USE_ENSEMBLE
    fig_dir = 'figures_ensemble_plots';
else
    fig_dir = sprintf('figures_snapshot_%d', SINGLE_SNAPSHOT);
end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Derived metrics 

if ndims(t_res_map) == 3
    tres_slice = squeeze(mean(t_res_map, 3, 'omitnan'));
else
    tres_slice = t_res_map;
end

% Wall sticking probability per (St, v0) cell
if ndims(outcome_map) == 3
    p_stick = squeeze(mean(double(outcome_map == 1), 3, 'omitnan'));
else
    p_stick = double(outcome_map == 1);
end

% Mean bounce count — NaN entries (never-stuck particles) → 0 before mean
bmap_clean = bounce_at_stick_map;
bmap_clean(isnan(bmap_clean)) = 0;
mean_bounces = squeeze(mean_omitnan_dim3(bmap_clean));

% Midplane crossing probability: fraction of (v0, dim3) cases with y_max > 0
midplane_prob = zeros(Ntau, 1);
for ia = 1:Ntau
    vals = ymax_slice(ia, :);          % use ensemble mean (or single snapshot slice)
    midplane_prob(ia) = mean(vals(:) > 0.0, 'omitnan');
end

% Mean penetration height vs St (averaged over all v0)
mean_ymax_vs_St = zeros(Ntau, 1);
for ia = 1:Ntau
    vals = ymax_slice(ia, :);
    mean_ymax_vs_St(ia) = mean(vals(:), 'omitnan');
end

%%  Figure 1: Maximum penetration height
fh = figure('Name','Maximum penetration height','Position',fig_size);
plot_log_heatmap(v0_norm, St_list, ymax_slice, 6);
colormap(parula);
cb = colorbar;
cb.Label.String   = 'y_{max}  (outer units)';
cb.Label.FontSize = 12;
clim([-1 1]);
xlabel('v_0 / U_{ref}',  'FontSize', 14);
ylabel(ref.St_label,     'FontSize', 14);
title(['Maximum penetration height  —  ' mode_label], 'FontSize', 13);
set(gca, 'Layer', 'top');
saveas(fh, fullfile(fig_dir, 'fig_01_ymax_map.png'));

%% Figure 2: Residence time
fh = figure('Name','Residence time','Position',fig_size);
plot_log_heatmap(v0_norm, St_list, tres_slice, 6);
colormap(parula);
cb = colorbar;
cb.Label.String   = 't_{res}  (outer units)';
cb.Label.FontSize = 12;
xlabel('v_0 / U_{ref}', 'FontSize', 14);
ylabel(ref.St_label,    'FontSize', 14);
title(['Residence time  —  ' mode_label], 'FontSize', 13);
set(gca, 'Layer', 'top');
saveas(fh, fullfile(fig_dir, 'fig_02_residence_time_map.png'));

%% Figure 3: Midplane crossing probability

log10St    = log10(St_list(:));
P_data     = midplane_prob(:);

[~, idx50] = min(abs(P_data - 0.5));
St50_init  = St_list(max(1, idx50));
k_init     = 5;
fit_ok     = false;

try
    fit_opts = optimset('Display', 'off');
    logistic = @(b, x) 1 ./ (1 + exp(-b(1) .* (x - log10(b(2)))));
    b_fit    = fminsearch(@(b) sum((logistic(b, log10St) - P_data).^2), ...
                          [k_init, St50_init], fit_opts);
    k_fit    = b_fit(1);
    St50_fit = b_fit(2);
    St_fine  = logspace(log10(St_list(1)), log10(St_list(end)), 300);
    P_fit    = logistic(b_fit, log10(St_fine));
    fit_ok   = true;
catch
    fprintf('Warning: logistic fit did not converge for Fig 3.\n');
end

fh = figure('Name','Midplane crossing probability','Position',fig_size);
ax3 = axes(fh);

if fit_ok
    fill(ax3, [St_fine, fliplr(St_fine)], [P_fit, zeros(1,numel(P_fit))], ...
         [0.85 0.2 0.1], 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    hold(ax3, 'on');
    semilogx(ax3, St_fine, P_fit, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 2.2, ...
             'DisplayName', sprintf('Logistic fit   St_{50} = %.3f,   k = %.1f', ...
             St50_fit, k_fit));
    hold(ax3, 'on');
end

semilogx(ax3, St_list, midplane_prob, 'o', 'MarkerSize', 7, ...
         'MarkerFaceColor', [0.15 0.4 0.75], 'MarkerEdgeColor', 'w', ...
         'LineWidth', 0.8, 'DisplayName', 'Simulation data');

if fit_ok
    xline(ax3, St50_fit, '--', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.2, ...
          'Alpha', 0.7, 'HandleVisibility', 'off');
    text(ax3, St50_fit * 1.08, 0.08, sprintf('St_{50} = %.3f', St50_fit), ...
         'FontSize', 10, 'Color', [0.85 0.2 0.1], 'FontName', 'Arial');
end

set(ax3, 'XScale', 'log');
xlabel(ax3, ref.St_label,        'FontSize', 14);
ylabel(ax3, 'P(y_{max} > 0)',    'FontSize', 14);
ylim(ax3,  [0 1]);
xlim(ax3,  [10^(log10(St_list(1)) - 0.1), 10^(log10(St_list(end)) + 0.1)]);
grid(ax3,  'on'); box(ax3, 'on');
legend(ax3, 'Interpreter', 'tex', 'Location', 'northwest', 'FontSize', 11);
title(ax3, ['Midplane crossing probability  —  ' mode_label], 'FontSize', 13);
saveas(fh, fullfile(fig_dir, 'fig_03_midplane_crossing_probability.png'));

%% Figure 4: Mean penetration scaling 

fh = figure('Name','Mean penetration scaling','Position',fig_size);
ax4 = axes(fh);

% Shaded midplane band
St_vec  = [St_list(:); flipud(St_list(:))];
band_hi = 0.05 * ones(Ntau, 1);
band_lo = -0.05 * ones(Ntau, 1);
fill(ax4, St_vec, [band_hi; flipud(band_lo)], [0.5 0.5 0.5], ...
     'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
set(ax4, 'XScale', 'log');
hold(ax4, 'on');

% Midplane reference
yline(ax4, 0, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
      'HandleVisibility', 'off');

% Data
semilogx(ax4, St_list, mean_ymax_vs_St ./ geom.H, 'o-', ...
         'MarkerSize', 6, 'MarkerFaceColor', [0.15 0.4 0.75], ...
         'MarkerEdgeColor', 'w', 'LineWidth', 1.8, 'Color', [0.15 0.4 0.75], ...
         'DisplayName', '\langle y_{max} \rangle / H');

xlabel(ax4, ref.St_label,               'FontSize', 14);
ylabel(ax4, '\langle y_{max} \rangle / H', 'FontSize', 14);
grid(ax4,  'on'); box(ax4, 'on');
xlim(ax4,  [10^(log10(St_list(1)) - 0.1), 10^(log10(St_list(end)) + 0.1)]);
title(ax4, ['Mean penetration height scaling  —  ' mode_label], 'FontSize', 13);
saveas(fh, fullfile(fig_dir, 'fig_04_mean_penetration_scaling.png'));

%% Figure 5: Regime map

regime_map_raw = zeros(size(ymax_slice));
regime_map_raw(ymax_slice > -0.5) = 1;
regime_map_raw(ymax_slice >  0.0) = 2;

% shadowing the ensemble-loaded alpha_partial / alpha_mid.
fig_fit_partial_ok = false;
fig_fit_mid_ok     = false;

v0_thresh_partial = NaN(Ntau, 1);
v0_thresh_mid     = NaN(Ntau, 1);
for ia = 1:Ntau
    row   = ymax_slice(ia, :);
    idx_p = find(row > -0.5, 1, 'first');
    idx_m = find(row >  0.0, 1, 'first');
    if ~isempty(idx_p), v0_thresh_partial(ia) = v0_norm(idx_p); end
    if ~isempty(idx_m), v0_thresh_mid(ia)     = v0_norm(idx_m); end
end

valid_p = ~isnan(v0_thresh_partial) & v0_thresh_partial > 0 & St_list(:) > 0;
if sum(valid_p) >= 4
    pp = polyfit(log(St_list(valid_p)), log(v0_thresh_partial(valid_p)), 1);
    fig_alpha_p = pp(1);
    fig_C_p     = exp(pp(2));
    fig_fit_partial_ok = true;
end

valid_m = ~isnan(v0_thresh_mid) & v0_thresh_mid > 0 & St_list(:) > 0;
if sum(valid_m) >= 4
    pm = polyfit(log(St_list(valid_m)), log(v0_thresh_mid(valid_m)), 1);
    fig_alpha_m = pm(1);
    fig_C_m     = exp(pm(2));
    fig_fit_mid_ok = true;
end

St_fit_line = logspace(log10(St_list(1)), log10(St_list(end)), 300);
log10St     = log10(St_list);

fh  = figure('Name','Regime map','Position',fig_size);
ax5 = axes(fh);

[v0f, log10Stf, regime_smooth] = smooth_grid(v0_norm, log10St, ...
                                             double(regime_map_raw), 6);
[Vf_grid, Sf_grid] = meshgrid(v0f, log10Stf);
pcolor(ax5, Vf_grid, Sf_grid, regime_smooth);
shading(ax5, 'flat');

colormap(ax5, [0.10 0.14 0.34;   % near-wall  (dark navy)
               0.26 0.55 0.76;   % partial    (steel blue)
               0.95 0.75 0.10]); % midplane   (amber)
hold(ax5, 'on');

% Raw-data contours — data-faithful boundary edges, not in legend
contour(ax5, v0_norm, log10St, ymax_slice, [-0.5 -0.5], ...
        'Color', 'w', 'LineWidth', 2.0, 'HandleVisibility', 'off');
contour(ax5, v0_norm, log10St, ymax_slice, [ 0.0  0.0], ...
        'Color', [0.15 0.15 0.15], 'LineStyle', '--', 'LineWidth', 1.8, ...
        'HandleVisibility', 'off');

% Power-law fits — only these appear in the legend
if fig_fit_partial_ok
    v0_fit_p = fig_C_p .* St_fit_line.^fig_alpha_p;
    plot(ax5, v0_fit_p, log10(St_fit_line), ':', 'Color', 'w', 'LineWidth', 1.8, ...
         'DisplayName', sprintf('v_0 \\propto St^{%.2f}  (partial)', fig_alpha_p));
end
if fig_fit_mid_ok
    v0_fit_m = fig_C_m .* St_fit_line.^fig_alpha_m;
    plot(ax5, v0_fit_m, log10(St_fit_line), ':', ...
         'Color', [0.15 0.15 0.15], 'LineWidth', 1.8, ...
         'DisplayName', sprintf('v_0 \\propto St^{%.2f}  (midplane)', fig_alpha_m));
end

% Y-axis: power-of-ten labels
decades   = floor(min(log10St)):ceil(max(log10St));
tick_vals = [];
for d = decades(1):decades(end)-1
    tick_vals = [tick_vals, d, d+0.5]; %#ok<AGROW>
end
tick_vals = [tick_vals, decades(end)];
tick_vals = tick_vals(tick_vals >= min(log10Stf) & tick_vals <= max(log10Stf));
ax5.YTick      = tick_vals;
ax5.YTickLabel = arrayfun(@(v) sprintf('10^{%g}', v), tick_vals, ...
                           'UniformOutput', false);

cb            = colorbar(ax5);
cb.Ticks      = [1/3, 1, 5/3];
cb.TickLabels = {'Near-wall', 'Partial penetration', 'Midplane crossing'};
cb.FontSize   = 11;
cb.TickLength = 0;
clim(ax5, [0 2]);

xlabel(ax5, 'v_0 / U_{ref}',                              'FontSize', 13);
ylabel(ax5, ref.St_label,                                  'FontSize', 13);
title(ax5,  ['Transport regime map  —  ' mode_label],     'FontSize', 13, ...
      'FontWeight', 'normal');
xlim(ax5, [v0_norm(1), v0_norm(end)]);
ylim(ax5, [min(log10Stf), max(log10Stf)]);
box(ax5, 'on');
ax5.Layer      = 'top';
ax5.FontSize   = 12;
ax5.FontName   = 'Arial';
ax5.XMinorTick = 'on';

if fig_fit_partial_ok || fig_fit_mid_ok
    leg5 = legend(ax5, 'Interpreter', 'tex', 'Location', 'northeast', 'FontSize', 10);
    leg5.Box       = 'off';
    leg5.TextColor = 'w';
end
saveas(fh, fullfile(fig_dir, 'fig_05_regime_map.png'));

% Console output
fprintf('\n--- Regime boundary fits ---\n');
if fig_fit_partial_ok
    fprintf('Partial penetration boundary: v0_crit = %.3f * St^(%.3f)\n', ...
            fig_C_p, fig_alpha_p);
end
if fig_fit_mid_ok
    fprintf('Midplane crossing boundary:   v0_crit = %.3f * St^(%.3f)\n', ...
            fig_C_m, fig_alpha_m);
end
fprintf('----------------------------\n\n');

%% Figure 6: Wall sticking probability 
fh = figure('Name','Wall sticking probability','Position',fig_size);
plot_log_heatmap(v0_norm, St_list, p_stick, 6);
colormap(parula);
clim([0 1]);
cb = colorbar;
cb.Label.String   = 'P(stick)';
cb.Label.FontSize = 12;
cb.FontSize       = 11;
xlabel('v_0 / U_{ref}', 'FontSize', 14);
ylabel(ref.St_label,    'FontSize', 14);
title(['Wall sticking probability  —  ' mode_label], 'FontSize', 13);
set(gca, 'Layer', 'top');
saveas(fh, fullfile(fig_dir, 'fig_06_stick_probability_map.png'));

%% Figure 7: Mean bounce count 
fh = figure('Name','Mean bounce count','Position',fig_size);
plot_log_heatmap(v0_norm, St_list, mean_bounces, 6);
colormap(parula);
cb = colorbar;
cb.Label.String   = 'Mean bounce count';
cb.Label.FontSize = 12;
cb.FontSize       = 11;
xlabel('v_0 / U_{ref}', 'FontSize', 14);
ylabel(ref.St_label,    'FontSize', 14);
title(['Mean bounce count at wall sticking  —  ' mode_label], 'FontSize', 13);
set(gca, 'Layer', 'top');
saveas(fh, fullfile(fig_dir, 'fig_07_mean_bounce_count_map.png'));

%% Figures 8 & 9: Trajectories 
if ~exist('has_trajectories','var') || ~has_trajectories
    fprintf('No trajectory data available — skipping Figs 8 and 9.\n');
else
    colors = lines(numel(trajectories));
    line_w = [2.2, 1.8, 1.4];   % thicker for low-St (more detail near wall)

    % Fig 8: x-y
    fh = figure('Name','Representative trajectories x-y','Position',fig_size);
    ax8 = axes(fh);
    set(ax8, 'Color', [0.97 0.97 0.97]);
    hold(ax8, 'on');
    for k = 1:numel(trajectories)
        s_hist = trajectories(k).s;
        lw = line_w(min(k, numel(line_w)));
        plot(ax8, s_hist(:,1), s_hist(:,2), '-', 'Color', colors(k,:), ...
             'LineWidth', lw, ...
             'DisplayName', sprintf('%s   \\tau_p = %.2e,  v_0 = %.2f', ...
             trajectories(k).label, trajectories(k).tau_p, trajectories(k).v0));
    end
    yline(ax8, geom.y_bottom, 'k-',  'LineWidth', 2.0, 'DisplayName', 'Wall');
    yline(ax8, geom.y_top,    'k-',  'LineWidth', 2.0, 'HandleVisibility', 'off');
    yline(ax8, 0,             '--',  'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, ...
          'DisplayName', 'Midplane');
    xlabel(ax8, 'x (outer units)', 'FontSize', 14);
    ylabel(ax8, 'y (outer units)', 'FontSize', 14);
    ylim(ax8, [geom.y_bottom - 0.06, geom.y_top + 0.06]);
    title(ax8, 'Representative particle trajectories', 'FontSize', 13);
    legend(ax8, 'Interpreter', 'none', 'Location', 'best', 'FontSize', 10);
    grid(ax8, 'on'); box(ax8, 'on');
    saveas(fh, fullfile(fig_dir, 'fig_08_representative_trajectories_xy.png'));

    % Fig 9: 3D
    fh = figure('Name','Representative trajectories 3D','Position',fig_size);
    ax9 = axes(fh);
    hold(ax9, 'on');
    for k = 1:numel(trajectories)
        s_hist = trajectories(k).s;
        lw = line_w(min(k, numel(line_w)));
        plot3(ax9, s_hist(:,1), s_hist(:,3), s_hist(:,2), '-', ...
              'Color', colors(k,:), 'LineWidth', lw, ...
              'DisplayName', trajectories(k).label);
    end
    xlabel(ax9, 'x', 'FontSize', 14);
    ylabel(ax9, 'z', 'FontSize', 14);
    zlabel(ax9, 'y', 'FontSize', 14);
    zlim(ax9,  [geom.y_bottom - 0.06, geom.y_top + 0.06]);
    title(ax9, 'Representative particle trajectories (3D)', 'FontSize', 13);
    legend(ax9, 'Interpreter', 'none', 'Location', 'best', 'FontSize', 10);
    grid(ax9, 'on'); box(ax9, 'on');
    view(ax9, 3);
    saveas(fh, fullfile(fig_dir, 'fig_09_representative_trajectories_3d.png'));
end

%% Figure 10: Deposition fraction vs St 

deposition_frac = zeros(Ntau, 1);
escape_frac     = zeros(Ntau, 1);
for ia = 1:Ntau
    outcomes = squeeze(outcome_map(ia,:,:));
    deposition_frac(ia) = mean(outcomes(:) == 1, 'omitnan');
    escape_frac(ia)     = mean(outcomes(:) == 0, 'omitnan');
end

fh   = figure('Name','Deposition and escape fraction vs St','Position',fig_size);
ax10 = axes(fh);

% Sum-to-1 reference
semilogx(ax10, St_list, deposition_frac + escape_frac, '-', ...
         'Color', [0.75 0.75 0.75], 'LineWidth', 0.9, ...
         'HandleVisibility', 'off');
hold(ax10, 'on');

semilogx(ax10, St_list, deposition_frac, 'o-', ...
         'Color', [0.85 0.2 0.1], 'MarkerSize', 6, ...
         'MarkerFaceColor', [0.85 0.2 0.1], 'MarkerEdgeColor', 'w', ...
         'LineWidth', 1.8, 'DisplayName', 'Wall deposition');

semilogx(ax10, St_list, escape_frac, 's-', ...
         'Color', [0.15 0.4 0.75], 'MarkerSize', 6, ...
         'MarkerFaceColor', [0.15 0.4 0.75], 'MarkerEdgeColor', 'w', ...
         'LineWidth', 1.8, 'DisplayName', 'Domain exit (non-deposition)');

xlabel(ax10, ref.St_label,                    'FontSize', 14);
ylabel(ax10, 'Fraction of injected particles', 'FontSize', 14);
ylim(ax10,  [0 1.05]);
xlim(ax10,  [10^(log10(St_list(1)) - 0.1), 10^(log10(St_list(end)) + 0.1)]);
grid(ax10, 'on'); box(ax10, 'on');
legend(ax10, 'Location', 'east', 'FontSize', 11);
title(ax10, ['Deposition and escape fraction  —  ' mode_label], 'FontSize', 13);
saveas(fh, fullfile(fig_dir, 'fig_10_deposition_fraction_vs_St.png'));

fprintf('Saved figures to %s\n', fig_dir);

%% Ensemble-only: variance map 
if USE_ENSEMBLE && ~isempty(ymax_std_slice)
    fh = figure('Name','Ensemble variance','Position',fig_size);
    plot_log_heatmap(v0_norm, St_list, ymax_std_slice, 6);
    colormap(plasma_colormap());
    cb = colorbar;
    cb.Label.String   = '\sigma(y_{max})  (outer units)';
    cb.Label.FontSize = 12;
    cb.FontSize       = 11;
    xlabel('v_0 / U_{ref}', 'FontSize', 14);
    ylabel(ref.St_label,    'FontSize', 14);
    title(['Ensemble std of y_{max}  —  11 snapshots'], 'FontSize', 13);
    set(gca, 'Layer', 'top');
    saveas(fh, fullfile(fig_dir, 'fig_11_ensemble_variance_map.png'));
end

%% Ensemble-only: convergence figure 
if USE_ENSEMBLE && ~isempty(running_mean_alpha_mid)
    fh   = figure('Name','Ensemble convergence','Position',[0 0 820 480]);
    ax_c = axes(fh);
    n_ax = 1:N_SNAPSHOTS;

    plot(ax_c, n_ax, running_mean_alpha_mid, 'o-', ...
         'Color', [0.85 0.2 0.1], 'MarkerFaceColor', [0.85 0.2 0.1], ...
         'MarkerEdgeColor', 'w', 'MarkerSize', 7, 'LineWidth', 1.8, ...
         'DisplayName', 'Running mean  \alpha_{mid}');
    hold(ax_c, 'on');

    plot(ax_c, n_ax, running_mean_alpha_partial, 's--', ...
         'Color', [0.15 0.4 0.75], 'MarkerFaceColor', [0.15 0.4 0.75], ...
         'MarkerEdgeColor', 'w', 'MarkerSize', 7, 'LineWidth', 1.8, ...
         'DisplayName', 'Running mean  \alpha_{partial}');

    % ±5% convergence band around final value
    final_val = running_mean_alpha_mid(end);
    patch(ax_c, [n_ax, fliplr(n_ax)], ...
          [final_val*1.05*ones(1,N_SNAPSHOTS), final_val*0.95*ones(1,N_SNAPSHOTS)], ...
          [0.85 0.2 0.1], 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');

    if ~isempty(converged_at)
        xline(ax_c, converged_at, 'k--', 'LineWidth', 1.4, ...
              'DisplayName', sprintf('Converged at n = %d', converged_at));
    end

    xlabel(ax_c, 'Number of snapshots in ensemble', 'FontSize', 14);
    ylabel(ax_c, '\alpha  (fitted boundary exponent)', 'FontSize', 14);
    xlim(ax_c,   [1 N_SNAPSHOTS]);
    xticks(ax_c, 1:N_SNAPSHOTS);
    grid(ax_c, 'on'); box(ax_c, 'on');
    legend(ax_c, 'Location', 'best', 'FontSize', 11);
    title(ax_c, 'Ensemble convergence: running mean of boundary exponent', ...
          'FontSize', 13);
    saveas(fh, fullfile(fig_dir, 'fig_12_ensemble_convergence.png'));
end

end  % end main function


%% Local helper: boundary fit
function [alpha, C] = fit_boundary_local(St_list, v0_norm, ymax_map, threshold, min_pts)
Ntau      = numel(St_list);
v0_thresh = nan(Ntau, 1);
for ia = 1:Ntau
    row = ymax_map(ia, :);
    idx = find(row > threshold, 1, 'first');
    if ~isempty(idx), v0_thresh(ia) = v0_norm(idx); end
end
valid = ~isnan(v0_thresh) & v0_thresh > 0 & St_list(:) > 0;
if sum(valid) < min_pts, alpha = NaN; C = NaN; return; end
p     = polyfit(log(St_list(valid)), log(v0_thresh(valid)), 1);
alpha = p(1);
C     = exp(p(2));
end


%% Local helper: plasma colormap
function cmap = plasma_colormap()
% Approximates the matplotlib 'plasma' colormap (256 colours, dark-to-bright)
t    = linspace(0, 1, 256)';
R    = 0.050 + 0.788.*t + 0.160.*t.^2;
G    = 0.030 + 0.090.*t.^0.5 - 0.020.*t;
B    = 0.527 + 0.100.*t - 0.600.*t.^2;
cmap = max(0, min(1, [R, G, B]));
end


%% Plot helpers
function plot_smooth_heatmap(x, y, Z, nFactor, varargin) %#ok<DEFNU>
[x_fine, y_fine, Z_fine] = smooth_grid(x, y, Z, nFactor);
imagesc(x_fine, y_fine, Z_fine, varargin{:});
set(gca, 'YDir', 'normal');
end

function plot_log_heatmap(x, y, Z, nFactor, varargin)
log10y = log10(y);
[x_fine, log10y_fine, Z_fine] = smooth_grid(x, log10y, Z, nFactor);
imagesc(x_fine, log10y_fine, Z_fine, varargin{:});
set(gca, 'YDir', 'normal');
ax = gca;
decades   = floor(min(log10y)):ceil(max(log10y));
tick_vals = [];
for d = decades(1):decades(end)-1
    tick_vals = [tick_vals, d, d+0.5]; %#ok<AGROW>
end
tick_vals = [tick_vals, decades(end)];
tick_vals = tick_vals(tick_vals >= min(log10y_fine) & tick_vals <= max(log10y_fine));
ax.YTick      = tick_vals;
ax.YTickLabel = arrayfun(@(v) sprintf('10^{%g}', v), tick_vals, ...
                          'UniformOutput', false);
end

function [x_fine, y_fine, Z_fine] = smooth_grid(x, y, Z, nFactor)
[X, Y]  = meshgrid(x, y);
x_fine  = linspace(min(x), max(x), nFactor * numel(x));
y_fine  = linspace(min(y), max(y), nFactor * numel(y));
[Xf, Yf] = meshgrid(x_fine, y_fine);
Z_fine  = interp2(X, Y, Z, Xf, Yf, 'spline');
end

function out = mean_omitnan_dim3(A)
[n1, n2, ~] = size(A);
out = nan(n1, n2);
for i = 1:n1
    for j = 1:n2
        vals = squeeze(A(i,j,:));
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            out(i,j) = mean(vals);
        end
    end
end
end