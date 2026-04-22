function Final_Laminar_Flow_Model_Plot()
clearvars; clc;

set(0,'DefaultAxesFontSize',13);
set(0,'DefaultLineLineWidth',1.6);
set(0,'DefaultFigureColor','w');

rng('default');   % reproducibility

%% Channel geometry and mean flow
H    = 1.0e-2;   % channel height [m]
L    = 0.20;     % channel length [m]
Uavg = 1.0;      % average streamwise velocity [m/s]

%% Parameter sweep  — outer-unit Stokes number St = tau_p * Uavg / H
tau_p_list = logspace(-3, 1, 41) .* (H/Uavg);   % dimensional tau_p [s]
v0_list    = linspace(0.1, 6.0, 41) .* Uavg;     % dimensional v0 [m/s]
St_list    = tau_p_list .* (Uavg/H);             % dimensionless St
v0_norm    = v0_list ./ Uavg;                     % dimensionless v0/Uavg

St_label = 'St = \tau_p U_{avg} / H';

Ntau = numel(tau_p_list);
Nv0  = numel(v0_list);

%% Wall model  — matches turbulent sweep exactly
p_stick      = 0.20;
N_bounce_max = 50;

%% Target — disabled
target.enabled = false;
target.x       = 0;
target.y       = 0;
target.eps     = 0;

%% Solver settings
y0      = 1e-12;              % start just above bottom wall [m]
tmax    = (L/Uavg) * 50;     % time cap [s]
RelTol  = 1e-7;
AbsTol  = 1e-10;
MaxStep = 0.02 * (H/Uavg);

%% Pre-allocate outputs
y_max_map   = nan(Ntau, Nv0);
outcome_map = zeros(Ntau, Nv0);
t_res_map   = nan(Ntau, Nv0);

%% Parameter sweep
fprintf('Running laminar sweep: %d x %d = %d cases...\n', Ntau, Nv0, Ntau*Nv0);
tic;

for ia = 1:Ntau
    tau_p = tau_p_list(ia);
    for ib = 1:Nv0
        v0 = v0_list(ib);

        [outcome, ~, ~, t_res, ~, ~, y_max] = single_particle_run( ...
            tau_p, v0, 0, ...
            H, L, Uavg, y0, ...
            p_stick, N_bounce_max, target, ...
            tmax, RelTol, AbsTol, MaxStep);

        outcome_map(ia,ib) = outcome;
        t_res_map(ia,ib)   = t_res;
        y_max_map(ia,ib)   = y_max;
    end
end

runtime = toc;
fprintf('Sweep complete in %.1f s.\n\n', runtime);

%% Normalize y_max to outer units matching turbulent script
ymax_norm = 2 .* (y_max_map ./ H) - 1;

%% Output folder
fig_dir = 'figures_laminar';
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
fig_size = [0 0 800 560];

%% FIGURE L1: y_max heatmap 
fh = figure('Name','L1 — Maximum penetration height','Position',fig_size);
plot_log_heatmap(v0_norm, St_list, ymax_norm, 6);
colormap(parula);
cb = colorbar;
cb.Label.String = 'y_{max} (outer units)';
clim([-1 1]);
xlabel('v_0 / U_{avg}');
ylabel(St_label);
title('Laminar: maximum penetration height y_{max}');
saveas(fh, fullfile(fig_dir,'fig_L1_ymax_map.png'));

%% FIGURE L2: Regime map with fitted power-law boundaries

regime_map = zeros(Ntau, Nv0);
regime_map(ymax_norm > -0.5) = 1;
regime_map(ymax_norm >  0.0) = 2;

% Boundary extraction
v0_thresh_partial = nan(Ntau,1);
v0_thresh_mid     = nan(Ntau,1);
for ia = 1:Ntau
    row   = ymax_norm(ia,:);
    idx_p = find(row > -0.5, 1, 'first');
    idx_m = find(row >  0.0, 1, 'first');
    if ~isempty(idx_p), v0_thresh_partial(ia) = v0_norm(idx_p); end
    if ~isempty(idx_m), v0_thresh_mid(ia)     = v0_norm(idx_m); end
end

% Power-law fits: v0_crit = C * St^alpha
fit_partial_ok = false;
fit_mid_ok     = false;

valid_p = ~isnan(v0_thresh_partial) & v0_thresh_partial > 0 & St_list(:) > 0;
if sum(valid_p) >= 4
    pp            = polyfit(log(St_list(valid_p)), log(v0_thresh_partial(valid_p)), 1);
    alpha_partial = pp(1);
    C_partial     = exp(pp(2));
    fit_partial_ok = true;
end

valid_m = ~isnan(v0_thresh_mid) & v0_thresh_mid > 0 & St_list(:) > 0;
if sum(valid_m) >= 4
    pm        = polyfit(log(St_list(valid_m)), log(v0_thresh_mid(valid_m)), 1);
    alpha_mid = pm(1);
    C_mid     = exp(pm(2));
    fit_mid_ok = true;
end

St_fit = logspace(log10(St_list(1)), log10(St_list(end)), 300);

fh  = figure('Name','L2 — Regime map','Position',fig_size);
ax  = axes(fh);
log10St = log10(St_list);
[v0f, log10Stf, regime_smooth] = smooth_grid(v0_norm, log10St, double(regime_map), 6);
[Vf, Sf] = meshgrid(v0f, log10Stf);
pcolor(ax, Vf, Sf, regime_smooth); shading(ax,'flat');
colormap(ax, [0.10 0.14 0.34;   % near-wall  (dark navy)
              0.26 0.55 0.76;   % partial    (steel blue)
              0.95 0.75 0.10]); % midplane   (amber)
hold(ax,'on');

% Raw-data contours
contour(ax, v0_norm, log10St, ymax_norm, [-0.5 -0.5], ...
        'Color','w', 'LineWidth', 2.0, 'HandleVisibility','off');
contour(ax, v0_norm, log10St, ymax_norm, [ 0.0  0.0], ...
        'Color',[0.15 0.15 0.15], 'LineStyle','--', 'LineWidth', 1.8, ...
        'HandleVisibility','off');

% Power-law fits 
if fit_partial_ok
    v0_fit_p = C_partial .* St_fit .^ alpha_partial;
    plot(ax, v0_fit_p, log10(St_fit), ':', 'Color','w', 'LineWidth', 1.8, ...
         'DisplayName', sprintf('v_0 \\propto St^{%.2f}  (partial)', alpha_partial));
end
if fit_mid_ok
    v0_fit_m = C_mid .* St_fit .^ alpha_mid;
    plot(ax, v0_fit_m, log10(St_fit), ':', 'Color',[0.15 0.15 0.15], 'LineWidth', 1.8, ...
         'DisplayName', sprintf('v_0 \\propto St^{%.2f}  (midplane)', alpha_mid));
end

% Y-axis: power-of-ten labels
decades   = floor(min(log10St)):ceil(max(log10St));
tick_vals = [];
for d = decades(1):decades(end)-1
    tick_vals = [tick_vals, d, d+0.5]; %#ok<AGROW>
end
tick_vals = [tick_vals, decades(end)];
tick_vals = tick_vals(tick_vals >= min(log10Stf) & tick_vals <= max(log10Stf));
ax.YTick      = tick_vals;
ax.YTickLabel = arrayfun(@(v) sprintf('10^{%g}',v), tick_vals, 'UniformOutput', false);

cb            = colorbar(ax);
cb.Ticks      = [1/3, 1, 5/3];
cb.TickLabels = {'Near-wall', 'Partial penetration', 'Midplane crossing'};
cb.FontSize   = 11;
cb.TickLength = 0;
clim(ax,[0 2]);
xlabel(ax, 'v_0 / U_{avg}',                  'FontSize', 13);
ylabel(ax, 'St = \tau_p \, U_{avg} / H',     'FontSize', 13);
title(ax,  'Laminar: transport regime map',   'FontSize', 13, 'FontWeight','normal');
xlim(ax, [v0_norm(1), v0_norm(end)]);
ylim(ax, [min(log10Stf), max(log10Stf)]);
box(ax,'on');
ax.Layer      = 'top';
ax.FontSize   = 12;
ax.FontName   = 'Arial';
ax.XMinorTick = 'on';
if fit_partial_ok || fit_mid_ok
    leg = legend(ax, 'Interpreter','tex', 'Location','northeast', 'FontSize', 10);
    leg.Box       = 'off';
    leg.TextColor = 'w';
end
saveas(fh, fullfile(fig_dir,'fig_L2_regime_map.png'));

fprintf('Laminar boundary fits:\n');
if fit_partial_ok
    fprintf('  Partial penetration: v0_crit = %.4f * St^(%.4f)\n', C_partial, alpha_partial);
end
if fit_mid_ok
    fprintf('  Midplane crossing:   v0_crit = %.4f * St^(%.4f)\n', C_mid, alpha_mid);
end

%% FIGURE L3: Midplane crossing probability vs St 
midplane_prob = zeros(Ntau,1);
for ia = 1:Ntau
    midplane_prob(ia) = mean(ymax_norm(ia,:) > 0, 'omitnan');
end

% Logistic fit: P = 1 / (1 + exp(-k*(log10(St) - log10(St50))))
log10St_vec = log10(St_list(:));
P_data      = midplane_prob(:);
[~, idx50]  = min(abs(P_data - 0.5));
St50_init   = St_list(max(1,idx50));
k_init      = 5;

logistic = @(b,x) 1 ./ (1 + exp(-b(1).*(x - log10(b(2)))));
fit_ok   = false;
try
    b_fit  = fminsearch(@(b) sum((logistic(b,log10St_vec) - P_data).^2), ...
                         [k_init, St50_init], optimset('Display','off'));
    St_fine = logspace(log10(St_list(1)), log10(St_list(end)), 300);
    P_fit   = logistic(b_fit, log10(St_fine));
    k_fit   = b_fit(1);
    St50_fit = b_fit(2);
    fit_ok  = true;
catch
    fprintf('Warning: logistic fit did not converge for Fig L3.\n');
end

fh = figure('Name','L3 — Midplane crossing probability','Position',fig_size);
semilogx(St_list, midplane_prob, 'o', 'MarkerSize',7, ...
         'MarkerFaceColor',[0.2 0.5 0.8], 'MarkerEdgeColor','none', ...
         'DisplayName','Laminar data');
hold on;
if fit_ok
    semilogx(St_fine, P_fit, '-', 'Color',[0.85 0.2 0.1], 'LineWidth',2, ...
             'DisplayName', sprintf('Logistic fit  St_{50}=%.3f,  k=%.1f', St50_fit, k_fit));
end
xlabel(St_label);
ylabel('P(y_{max} > 0)');
ylim([0 1]);
xlim([10^(log10(St_list(1))-0.1), 10^(log10(St_list(end))+0.1)]);
grid on; box on;
legend('Location','northwest');
title('Laminar: midplane crossing probability vs Stokes number');
saveas(fh, fullfile(fig_dir,'fig_L3_midplane_probability.png'));

if fit_ok
    fprintf('Laminar logistic fit: St50 = %.4f,  k = %.2f\n', St50_fit, k_fit);
end

%% FIGURE L4: Analytical consistency check (R5 verification)

[TAU, V0] = meshgrid(tau_p_list, v0_list);   % [Nv0 x Ntau]
ymax_analytical_raw = V0 .* TAU;             % [Nv0 x Ntau] in metres
ymax_analytical_norm = 2.*(ymax_analytical_raw./H) - 1;
ymax_analytical_norm = ymax_analytical_norm.';  % transpose to [Ntau x Nv0]

ymax_analytical_unclamped = ymax_analytical_norm;
ymax_analytical_norm      = min(max(ymax_analytical_norm, -1), 1);

abs_err_norm = abs(ymax_norm - ymax_analytical_norm) * 100;   % % of H/2

valid_mask = (ymax_analytical_norm > -0.95) & ...
             (ymax_analytical_norm <  0.95) & ...
             isfinite(ymax_norm) & isfinite(ymax_analytical_norm) & ...
             (St_list(:) < 3.0);

rel_err_pct       = abs_err_norm;
rel_err_pct_valid = rel_err_pct(valid_mask & isfinite(rel_err_pct));

max_err  = max(rel_err_pct_valid);
mean_err = mean(rel_err_pct_valid);

combined_mask = valid_mask & isfinite(rel_err_pct);
err_full = rel_err_pct;
err_full(~combined_mask) = -Inf;          % hide invalid cells
[~, sort_idx] = sort(err_full(:), 'descend');
fprintf('\nTop 5 error cases (St, v0, ymax_num, ymax_ana_clamped, err%%):\n');
for kk = 1:5
    [ri, ci] = ind2sub(size(rel_err_pct), sort_idx(kk));
    fprintf('  St=%.4f  v0=%.2f  num=%.4f  ana_clamp=%.4f  ana_raw=%.4f  err=%.4f%%\n', ...
        St_list(ri), v0_norm(ci), ...
        ymax_norm(ri,ci), ymax_analytical_norm(ri,ci), ...
        ymax_analytical_unclamped(ri,ci), rel_err_pct(ri,ci));
end

fprintf('\nAnalytical consistency (R5):\n');
fprintf('  Valid cases (analytical prediction within channel): %d of %d\n', ...
    sum(valid_mask(:)), numel(valid_mask));
fprintf('  Max error (as %% of channel half-width):  %.2f%%\n', max_err);
fprintf('  Mean error (as %% of channel half-width): %.2f%%\n', mean_err);
if max_err < 2.0
    fprintf('  R5 SATISFIED (< 2%% criterion)\n\n');
else
    fprintf('  R5 NOT SATISFIED — check event detection\n\n');
end

rel_err_display = rel_err_pct;
rel_err_display(~valid_mask) = NaN;

clim_max = max(prctile(rel_err_pct_valid, 99), 0.1);   % floor at 0.1%

log10St_err  = log10(St_list);
log10Stf_err = log10St_err;
v0f_err      = v0_norm;

SENTINEL     = clim_max * 2;          % safely above colour range
err_display2 = rel_err_display;       % NaN where excluded
err_display2(~valid_mask) = SENTINEL; % replace NaN with sentinel


n_steps  = 240;
navy     = [0.10, 0.14, 0.34];
sblue    = [0.26, 0.55, 0.76];
amber    = [0.95, 0.75, 0.10];

% Two-segment interpolation
half     = floor(n_steps/2);
seg1     = [linspace(navy(1),sblue(1),half)', ...
            linspace(navy(2),sblue(2),half)', ...
            linspace(navy(3),sblue(3),half)'];
seg2     = [linspace(sblue(1),amber(1),n_steps-half)', ...
            linspace(sblue(2),amber(2),n_steps-half)', ...
            linspace(sblue(3),amber(3),n_steps-half)'];
cmap_valid = [seg1; seg2];

n_grey   = 20;
cmap_grey = repmat([0.62 0.62 0.62], n_grey, 1);   % mid-grey for excluded
cmap_full = [cmap_valid; cmap_grey];

% Draw 
fh  = figure('Name','L4 — Analytical consistency','Position',fig_size,'Color','w');
ax4 = axes(fh);
imagesc(ax4, v0_norm, log10St_err, err_display2);
set(ax4,'YDir','normal');
colormap(ax4, cmap_full);
clim(ax4, [0, SENTINEL]);
hold(ax4,'on');

% Boundary contour
St3_log = log10(3.0);
xline(ax4, NaN);   % dummy — use yline on the log-St axis instead
plot(ax4, [v0_norm(1), v0_norm(end)], [St3_log, St3_log], ...
     'w-', 'LineWidth', 1.2, 'HandleVisibility','off');

% Y-axis: power-of-ten labels 
decades_err = floor(min(log10St_err)):ceil(max(log10St_err));
tv_err = [];
for d = decades_err(1):decades_err(end)-1
    tv_err = [tv_err, d, d+0.5]; %#ok<AGROW>
end
tv_err = [tv_err, decades_err(end)];
tv_err = tv_err(tv_err >= min(log10Stf_err) & tv_err <= max(log10Stf_err));
ax4.YTick      = tv_err;
ax4.YTickLabel = arrayfun(@(v) sprintf('10^{%g}',v), tv_err, 'UniformOutput',false);

% Colourbar
cb4 = colorbar(ax4);
cb4.Label.String   = 'Relative error  (% of H/2)';
cb4.Label.FontSize = 11;
cb4.FontSize       = 11;
cb4.TickLength     = 0;
cb4.Limits         = [0, clim_max];
n_ticks            = 5;
cb4.Ticks          = linspace(0, clim_max, n_ticks);
cb4.TickLabels     = arrayfun(@(v) sprintf('%.3f%%', v), cb4.Ticks, ...
                               'UniformOutput', false);

% Region labels 
text(ax4, 0.97, 0.97, 'Outside R5 scope  (St \geq 3)', ...
     'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
     'FontSize', 10, 'Color','w', 'FontName','Arial');

% Max / mean annotation 
text(ax4, 0.97, 0.03, sprintf('Max = %.2f%%   Mean = %.2f%%', max_err, mean_err), ...
     'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
     'FontSize', 10, 'Color','w', 'FontName','Arial');

% Axes cosmetics 
xlabel(ax4, 'v_0 / U_{avg}',                       'FontSize', 13);
ylabel(ax4, 'St = \tau_p \, U_{avg} / H',          'FontSize', 13);
title(ax4,  'Laminar: analytical consistency (R5)', 'FontSize', 13, 'FontWeight','normal');
xlim(ax4, [v0_norm(1), v0_norm(end)]);
ylim(ax4, [min(log10Stf_err), max(log10Stf_err)]);
box(ax4,'on');
ax4.Layer      = 'top';
ax4.FontSize   = 12;
ax4.FontName   = 'Arial';
ax4.XMinorTick = 'on';
legend(ax4,'off');
saveas(fh, fullfile(fig_dir,'fig_L4_analytical_consistency.png'));

%% FIGURE L5: Representative trajectories x-y 
traj_cases = struct( ...
    'label', {'low\_St\_near\_wall', 'mid\_St\_transitional', 'high\_St\_ballistic'}, ...
    'tau_p', {tau_p_list(6), tau_p_list(21), tau_p_list(31)}, ...
    'v0',    {v0_list(ceil(Nv0*0.6)), v0_list(ceil(Nv0*0.6)), v0_list(ceil(Nv0*0.4))});

colors = lines(numel(traj_cases));
fh = figure('Name','L5 — Representative trajectories','Position',fig_size);
hold on;

for k = 1:numel(traj_cases)
    tau_p = traj_cases(k).tau_p;
    v0    = traj_cases(k).v0;
    St_k  = tau_p * Uavg / H;
    v0_k  = v0 / Uavg;

    [t_hist, s_hist] = single_particle_run_history( ...
        tau_p, v0, 0, ...
        H, L, Uavg, y0, ...
        p_stick, N_bounce_max, target, ...
        tmax, RelTol, AbsTol, MaxStep);

    % Normalize x and y to outer units
    x_norm = s_hist(:,1) ./ H;          % x/H
    y_norm = 2.*(s_hist(:,2)./H) - 1;   % [-1,+1]

    plot(x_norm, y_norm, '-', 'Color', colors(k,:), 'LineWidth', 1.6, ...
         'DisplayName', sprintf('St=%.3f,  v_0/U=%.2f', St_k, v0_k));
end

% Channel walls and midplane
yline(-1, 'k-',  'LineWidth', 1.2, 'DisplayName', 'Wall');
yline( 1, 'k-',  'LineWidth', 1.2, 'HandleVisibility', 'off');
yline( 0, 'k--', 'LineWidth', 0.9, 'DisplayName', 'Midplane');

xlabel('x / H');
ylabel('y  (outer units,  -1 = bottom wall,  +1 = top wall)');
ylim([-1.05, 1.05]);
title('Laminar: representative particle trajectories');
legend('Interpreter','none','Location','best','FontSize',10);
grid on; box on;
saveas(fh, fullfile(fig_dir,'fig_L5_representative_trajectories.png'));

fprintf('Saved all figures to %s/\n', fig_dir);

end  % end main function

%% Single particle run 
function [outcome, x_final, y_final, t_total, crossed_target, ...
          bounce_at_stick, y_max] = single_particle_run( ...
          tau_p, v0_perp, x0, ...
          H, L, Uavg, y0, ...
          p_stick, N_bounce_max, target, ...
          tmax, RelTol, AbsTol, MaxStep)

outcome         = 0;
x_final         = NaN;
y_final         = NaN;
t_total         = 0.0;
crossed_target  = false;
bounce_at_stick = NaN;
y_max           = y0;

s0      = [x0; y0; 0.0; v0_perp];
bounced = 0;

while true
    rhs     = @(t,s) particle_rhs(s, tau_p, H, Uavg);
    eventsf = @(t,s) particle_events(s, H, L, target);
    opts    = odeset('RelTol',RelTol,'AbsTol',AbsTol, ...
                     'MaxStep',MaxStep,'Events',eventsf);

    [t_sol, s_sol, ~, se, ie] = ode45(rhs, [0, tmax - t_total], s0, opts);

    t_total = t_total + t_sol(end);
    y_max   = max(y_max, max(s_sol(:,2)));

    if isempty(ie)
        outcome = 0;
        x_final = s_sol(end,1); y_final = s_sol(end,2);
        break;
    end

    term_mask = (ie == 1) | (ie == 2) | (ie == 3);

    if ~any(term_mask)
        if t_total >= tmax
            outcome = 0;
            x_final = s_sol(end,1); y_final = s_sol(end,2);
            break;
        else
            s0 = s_sol(end,:).';
            continue;
        end
    end

    idx_last = find(term_mask, 1, 'last');
    ev = ie(idx_last);
    sv = se(idx_last,:);

    if ev == 3
        outcome = 0;
        x_final = sv(1); y_final = sv(2);
        break;
    end

    if rand() <= p_stick
        outcome         = 1;
        x_final         = sv(1); y_final = sv(2);
        bounce_at_stick = bounced;
        break;
    end

    bounced = bounced + 1;
    if bounced >= N_bounce_max
        outcome         = 1;
        x_final         = sv(1); y_final = sv(2);
        bounce_at_stick = bounced;
        break;
    end

    xR = sv(1); yR = sv(2); vxR = sv(3); vyR = sv(4);
    if ev == 1
        vyR = abs(vyR);  yR = max(yR, 0)   + 1e-12;
    else
        vyR = -abs(vyR); yR = min(yR, H) - 1e-12;
    end
    s0 = [xR; yR; vxR; vyR];
end

end


%% Single particle run with full history (for trajectory plots)
function [t_hist, s_hist] = single_particle_run_history( ...
          tau_p, v0_perp, x0, ...
          H, L, Uavg, y0, ...
          p_stick, N_bounce_max, target, ...
          tmax, RelTol, AbsTol, MaxStep)

s0      = [x0; y0; 0.0; v0_perp];
bounced = 0;
t_total = 0.0;
t_segs  = {};
s_segs  = {};

while true
    rhs     = @(t,s) particle_rhs(s, tau_p, H, Uavg);
    eventsf = @(t,s) particle_events(s, H, L, target);
    opts    = odeset('RelTol',RelTol,'AbsTol',AbsTol, ...
                     'MaxStep',MaxStep,'Events',eventsf);

    [t_sol, s_sol, ~, se, ie] = ode45(rhs, [0, tmax - t_total], s0, opts);

    if isempty(t_segs)
        t_segs{end+1} = t_sol + t_total;           %#ok<AGROW>
        s_segs{end+1} = s_sol;                      %#ok<AGROW>
    else
        t_segs{end+1} = t_sol(2:end) + t_total;    %#ok<AGROW>
        s_segs{end+1} = s_sol(2:end,:);             %#ok<AGROW>
    end

    t_total = t_total + t_sol(end);

    if isempty(ie)
        break;
    end

    term_mask = (ie == 1) | (ie == 2) | (ie == 3);
    if ~any(term_mask)
        if t_total >= tmax, break; end
        s0 = s_sol(end,:).';
        continue;
    end

    idx_last = find(term_mask, 1, 'last');
    ev = ie(idx_last);
    sv = se(idx_last,:);

    if ev == 3, break; end

    if rand() <= p_stick, break; end

    bounced = bounced + 1;
    if bounced >= N_bounce_max, break; end

    xR = sv(1); yR = sv(2); vxR = sv(3); vyR = sv(4);
    if ev == 1
        vyR = abs(vyR);  yR = max(yR, 0)   + 1e-12;
    else
        vyR = -abs(vyR); yR = min(yR, H) - 1e-12;
    end
    s0 = [xR; yR; vxR; vyR];
end

t_hist = vertcat(t_segs{:});
s_hist = vertcat(s_segs{:});

end


%% ODE RHS: Stokes drag in Poiseuille flow (no gravity)
function ds = particle_rhs(s, tau_p, H, Uavg)

y  = s(2);
vx = s(3);
vy = s(4);

yclamp = min(max(y, 0), H);
u  = 6 * Uavg * (yclamp/H) * (1 - yclamp/H);
v  = 0;

ax = (u - vx) / tau_p;
ay = (v - vy) / tau_p;

ds = [vx; vy; ax; ay];

end


%% Events: bottom wall, top wall, outlet (target disabled)
function [value, isterminal, direction] = particle_events(s, H, L, target)

x = s(1); y = s(2);

v_bottom = y;
v_top    = H - y;
v_out    = L - x;

if target.enabled
    dx = x - target.x; dy = y - target.y;
    v_tgt = sqrt(dx*dx + dy*dy) - target.eps;
    is_term_target = 0;
else
    v_tgt          = 1;
    is_term_target = 0;
end

value      = [v_bottom; v_top; v_out; v_tgt];
isterminal = [1; 1; 1; is_term_target];
direction  = [-1; -1; -1; -1];

end


%% Plot helpers

function plot_log_heatmap(x, y, Z, nFactor)
log10y = log10(y);
[x_fine, log10y_fine, Z_fine] = smooth_grid(x, log10y, Z, nFactor);
imagesc(x_fine, log10y_fine, Z_fine);
set(gca,'YDir','normal');
ax = gca;
decades   = floor(min(log10y)):ceil(max(log10y));
tick_vals = [];
for d = decades(1):decades(end)-1
    tick_vals = [tick_vals, d, d+0.5]; %#ok<AGROW>
end
tick_vals = [tick_vals, decades(end)];
tick_vals = tick_vals(tick_vals >= min(log10y_fine) & tick_vals <= max(log10y_fine));
ax.YTick      = tick_vals;
ax.YTickLabel = arrayfun(@(v) sprintf('10^{%g}',v), tick_vals, 'UniformOutput',false);
end

function [x_fine, y_fine, Z_fine] = smooth_grid(x, y, Z, nFactor)
[X, Y]   = meshgrid(x, y);
x_fine   = linspace(min(x), max(x), nFactor*numel(x));
y_fine   = linspace(min(y), max(y), nFactor*numel(y));
[Xf, Yf] = meshgrid(x_fine, y_fine);
Z_fine   = interp2(X, Y, Z, Xf, Yf, 'spline');
end
