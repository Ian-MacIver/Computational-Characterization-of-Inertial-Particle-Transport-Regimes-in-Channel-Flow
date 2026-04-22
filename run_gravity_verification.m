function run_gravity_verification()

clearvars; clc; rng('default');

set(0,'DefaultAxesFontSize',13);
set(0,'DefaultLineLineWidth',1.6);
set(0,'DefaultFigureColor','w');

%% Channel geometry — MUST match Final_Laminar_Flow_Model_Plot.m exactly
H    = 1.0e-2;   % channel height [m]
L    = 0.20;     % channel length [m]
Uavg = 1.0;      % average streamwise velocity [m/s]

%% Physicial Constants
g_phys  = 9.81;                       % gravitational acceleration [m/s²]
Gamma   = g_phys * H / Uavg^2;        % dimensionless gravitational parameter
                                       % Gamma = g*H/Uref² — outer units
fprintf('=============================================================\n');
fprintf('  Gravity verification (R6)\n');
fprintf('=============================================================\n');
fprintf('Channel: H = %.4g m,  L = %.2f m,  Uavg = %.2f m/s\n', H, L, Uavg);
fprintf('Gamma = g*H/Uref^2 = %.4f  (dimensionless gravitational parameter)\n', Gamma);
fprintf('Gravity shift per unit St (outer units): 2g*H/Uref^2 * Uref/H = %.4f\n\n', 2*g_phys*H/Uavg^2);

%%  Solver settings — MUST match Final_Laminar_Flow_Model_Plot.m exactly
y0_init      = 1e-12;              % initial height just above bottom wall [m]
tmax         = (L/Uavg) * 50;     % time cap [s]
RelTol       = 1e-7;
AbsTol       = 1e-10;
MaxStep      = 0.02 * (H/Uavg);

%% Wall Model
p_stick      = 0.20;
N_bounce_max = 50;

%% Target Settings
target.enabled = false;
target.x       = 0;
target.y       = 0;
target.eps     = 0;



tau_p_cases = [H/Uavg * 0.050;   % St = 0.050
               H/Uavg * 0.100;   % St = 0.100
               H/Uavg * 0.150;   % St = 0.150
               H/Uavg * 0.200;   % St = 0.200
               H/Uavg * 0.250];  % St = 0.250

% v0 = 0.15 * Uavg / St → ymax_norm = -0.70 (deep partial penetration)
v0_cases = [0.15/0.050 * Uavg;   % St=0.050 → v0 = 3.00
            0.15/0.100 * Uavg;   % St=0.100 → v0 = 1.50
            0.15/0.150 * Uavg;   % St=0.150 → v0 = 1.00
            0.15/0.200 * Uavg;   % St=0.200 → v0 = 0.75
            0.15/0.250 * Uavg];  % St=0.250 → v0 = 0.60

N_cases = numel(tau_p_cases);

% Verify effective injection velocity is positive for all cases
fprintf('Pre-check: effective injection velocity v0_eff = v0 - g*tau_p > 0 ?\n');
all_ok = true;
for k = 1:N_cases
    v0_eff = v0_cases(k) - g_phys * tau_p_cases(k);
    St_k   = tau_p_cases(k) * Uavg / H;
    ok_str = 'OK';
    if v0_eff <= 0
        ok_str  = '*** NEGATIVE — particle never rises! Increase v0 ***';
        all_ok  = false;
    end
    fprintf('  Case %d: St=%.4f  v0=%.4f  g*tau_p=%.6f  v0_eff=%.6f  [%s]\n', ...
        k, St_k, v0_cases(k), g_phys*tau_p_cases(k), v0_eff, ok_str);
end
if ~all_ok
    error('One or more cases have non-positive effective injection velocity. Increase v0.');
end
fprintf('NOTE: Verification runs use uniform streamwise flow (u_f = Uavg = const)\n');
fprintf('      This eliminates Poiseuille x-y shear coupling so the analytical\n');
fprintf('      formula y_max = (v0 - g*tau_p)*tau_p holds to solver tolerance.\n');
fprintf('      The gravity translation argument is independent of the streamwise\n');
fprintf('      profile — it acts only on the wall-normal equation.\n\n'); 
fprintf('Boundary position check (cases must be BELOW midplane boundary):\n');
fprintf('  %-6s %-8s %-10s %-14s %-14s %-12s %-10s\n', ...
    'Case', 'St', 'v0/Uref', 'v0_partial', 'v0_midplane', 'ymax_pred', 'Status');
fprintf('  %s\n', repmat('-', 1, 74));
C_mid_chk     = 0.6408;   alpha_mid_chk     = -0.8500;
C_partial_chk = 0.3700;   alpha_partial_chk = -0.8089;
for k = 1:N_cases
    St_k      = tau_p_cases(k) * Uavg / H;
    v0_k      = v0_cases(k) / Uavg;
    v0_p      = C_partial_chk * St_k^alpha_partial_chk;
    v0_m      = C_mid_chk     * St_k^alpha_mid_chk;
    ymax_pred_k = 2*((v0_cases(k) - g_phys*tau_p_cases(k))*tau_p_cases(k)/H) - 1;
    if v0_k < v0_p
        status = 'near-wall';
    elseif v0_k < v0_m
        status = 'partial ✓';
    else
        status = '*** ABOVE MIDPLANE — will hit wall! ***';
    end
    fprintf('  %-6d %-8.4f %-10.4f %-14.4f %-14.4f %-12.4f %-10s\n', ...
        k, St_k, v0_k, v0_p, v0_m, ymax_pred_k, status);
end
fprintf('\n');

p_stick_verif  = 1.0;
N_bounce_verif = 0;
L_verif        = 1e6 * H;     % no outlet termination
tmax_verif     = 200 * max(tau_p_cases);  % >> 2*tau_p for all cases
fprintf('%-6s %-10s %-10s %-14s %-14s %-14s %-12s\n', ...
    'Case', 'St', 'v0/Uavg', ...
    'ymax,norm(g=0)', 'ymax,norm(g>0)', 'ymax,norm(pred)', 'Error (%)');
fprintf('%s\n', repmat('-', 1, 80));

max_error = 0;
% Store results: [St, v0_norm, ymax_off_norm, ymax_on_norm, ymax_pred_norm, err]
results   = zeros(N_cases, 6);

for k = 1:N_cases
    tau_p = tau_p_cases(k);
    v0    = v0_cases(k);
    St    = tau_p * Uavg / H;

    % Gravity Disabled
    [~, ~, ~, ~, ~, ~, ymax_off_dim] = single_particle_run_gravity( ...
        tau_p, v0, 0, ...
        H, L_verif, Uavg, y0_init, ...
        p_stick_verif, N_bounce_verif, target, ...
        tmax_verif, RelTol, AbsTol, MaxStep, ...
        false, g_phys, true);   % useUniformFlow=true: decouples x-y exactly

    % Gravity Enabled
    [~, ~, ~, ~, ~, ~, ymax_on_dim] = single_particle_run_gravity( ...
        tau_p, v0, 0, ...
        H, L_verif, Uavg, y0_init, ...
        p_stick_verif, N_bounce_verif, target, ...
        tmax_verif, RelTol, AbsTol, MaxStep, ...
        true, g_phys, true);    % useUniformFlow=true: same flow as gravity-off run

    % Prediction (Analytical)
    % y_max,predicted = (v0 - g*tau_p) * tau_p   [metres]
    ymax_pred_dim = (v0 - g_phys * tau_p) * tau_p;

    % Outer-unit normalisation [-1, +1] 
    % y_norm = 2*(y_dim/H) - 1
    ymax_off_norm  = 2 * (ymax_off_dim  / H) - 1;
    ymax_on_norm   = 2 * (ymax_on_dim   / H) - 1;
    ymax_pred_norm = 2 * (ymax_pred_dim / H) - 1;

    % Rel. Error
    if abs(ymax_pred_norm) > 1e-10
        err = abs(ymax_on_norm - ymax_pred_norm) / abs(ymax_pred_norm) * 100;
    else
        err = NaN;
        fprintf('  WARNING: Case %d analytical prediction near zero — error undefined\n', k);
    end

    max_error = max(max_error, err);

    results(k,:) = [St, v0/Uavg, ymax_off_norm, ymax_on_norm, ymax_pred_norm, err];

    fprintf('%-6d %-10.5f %-10.4f %-14.6f %-14.6f %-14.6f %-12.4f\n', ...
        k, St, v0/Uavg, ymax_off_norm, ymax_on_norm, ymax_pred_norm, err);
end

fprintf('%s\n', repmat('-', 1, 80));
fprintf('Maximum relative error across all cases: %.4f%%\n\n', max_error);

%% R6 Check
if max_error < 2.0
    fprintf('✓ R6 SATISFIED: all five cases agree with analytical prediction\n');
    fprintf('  within %.4f%% (< 2%% criterion)\n\n', max_error);
else
    fprintf('✗ R6 NOT SATISFIED: max error %.4f%% exceeds 2%% — check gravity RHS\n\n', max_error);
end

%% Table printout 
fprintf('=============================================================\n');
fprintf('  Gravity verification (outer-unit normalisation)\n');
fprintf('  y_norm = 2*(y/H) - 1  ∈ [-1, +1];  injection wall = -1, midplane = 0\n');
fprintf('=============================================================\n');
fprintf('%-6s %-10s %-10s %-18s %-18s %-18s %-12s\n', ...
    'Case', 'St', 'v0/Uref', ...
    'y_max (no gravity)', 'y_max (pred, grav)', 'y_max (num, grav)', 'Rel. error');
fprintf('%s\n', repmat('-', 1, 96));
for k = 1:N_cases
    fprintf('%-6d %-10.5f %-10.4f %-18.6f %-18.6f %-18.6f %-12s\n', ...
        k, results(k,1), results(k,2), ...
        results(k,3), results(k,5), results(k,4), ...
        sprintf('%.4f%%', results(k,6)));
end
fprintf('%s\n\n', repmat('-', 1, 96));

%% Results for Post-Processing
save('gravity_verification_results.mat', ...
    'results', 'tau_p_cases', 'v0_cases', 'H', 'Uavg', 'g_phys', 'Gamma', 'max_error');
fprintf('Results saved to gravity_verification_results.mat\n\n');

%% Laminar Boundary Fit
% Typical laminar values: alpha ~ -1.00 for both boundaries.
C_partial     = 0.30;     % ← replace with your laminar fit output
alpha_partial = -1.00;    % ← replace with your laminar fit output
C_mid         = 0.60;     % ← replace with your laminar fit output
alpha_mid     = -1.00;    % ← replace with your laminar fit output

Gamma_outer = Gamma;   % = g * H / Uref^2 in your normalization

fprintf('Using laminar boundary fits:\n');
fprintf('  Partial:  v0_crit = %.4f * St^(%.4f)\n', C_partial, alpha_partial);
fprintf('  Midplane: v0_crit = %.4f * St^(%.4f)\n\n', C_mid, alpha_mid);

%% FIGURE 1: Regime boundary translation 

St_plot = logspace(-2, 1, 500);

% Gravity-free boundaries
v0_partial_no_g = C_partial .* St_plot .^ alpha_partial;
v0_mid_no_g     = C_mid     .* St_plot .^ alpha_mid;

% Gravity-shifted boundaries
v0_partial_grav = v0_partial_no_g + Gamma_outer .* St_plot;
v0_mid_grav     = v0_mid_no_g     + Gamma_outer .* St_plot;

fig_size = [0 0 860 560];
v0_lo = 0.08;  v0_hi = 12;   % y-axis range

fh1 = figure('Name', 'Fig 4.12 — Gravity regime boundary translation', ...
             'Position', fig_size, 'Color', 'w');
ax1 = axes(fh1);
set(ax1, 'XScale', 'log', 'YScale', 'log'); hold(ax1, 'on');

% Background regime zones (no gravity) 
St_fill = [St_plot, fliplr(St_plot)];

% Near-wall zone (below partial penetration line)
fill(ax1, St_fill, [max(v0_partial_no_g, v0_lo), fliplr(repmat(v0_lo, size(St_plot)))], ...
     [0.10 0.14 0.34], 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Partial penetration zone (between partial and midplane)
fill(ax1, St_fill, [v0_mid_no_g, fliplr(v0_partial_no_g)], ...
     [0.26 0.55 0.76], 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Midplane crossing zone (above midplane line, up to axis top)
fill(ax1, St_fill, [repmat(v0_hi, size(St_plot)), fliplr(v0_mid_no_g)], ...
     [0.95 0.75 0.10], 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Gravity-free boundaries 
h_p_ng = loglog(ax1, St_plot, v0_partial_no_g, '-', ...
    'Color', [0.15 0.35 0.75], 'LineWidth', 2.2, ...
    'DisplayName', sprintf('Partial boundary, g = 0  (v_0 \\propto St^{%.2f})', alpha_partial));
h_m_ng = loglog(ax1, St_plot, v0_mid_no_g, '--', ...
    'Color', [0.15 0.35 0.75], 'LineWidth', 2.2, ...
    'DisplayName', sprintf('Midplane boundary, g = 0  (v_0 \\propto St^{%.2f})', alpha_mid));

% Gravity-shifted boundaries 
h_p_g = loglog(ax1, St_plot, v0_partial_grav, '-', ...
    'Color', [0.80 0.15 0.10], 'LineWidth', 2.2, ...
    'DisplayName', sprintf('Partial boundary, g > 0  (shift = \\Gamma\\cdotSt,  \\Gamma = %.4f)', Gamma_outer));
h_m_g = loglog(ax1, St_plot, v0_mid_grav, '--', ...
    'Color', [0.80 0.15 0.10], 'LineWidth', 2.2, ...
    'DisplayName', 'Midplane boundary, g > 0');

% Shade the translation gap for the midplane boundary 
v0_top_clip  = min(v0_mid_grav, v0_hi);
v0_bot_clip  = max(v0_mid_no_g, v0_lo);
fill(ax1, St_fill, [v0_top_clip, fliplr(v0_bot_clip)], ...
     [0.90 0.30 0.20], 'FaceAlpha', 0.13, 'EdgeColor', 'none', 'HandleVisibility', 'off');

St_arrow = 2.0;
v0_arrow_bot = C_mid * St_arrow^alpha_mid;
v0_arrow_top = v0_arrow_bot + Gamma_outer * St_arrow;
if v0_arrow_top <= v0_hi && v0_arrow_bot >= v0_lo
    annotation(fh1, 'doublearrow', ...
        [0 0], [0 0], 'Color', [0.80 0.15 0.10]);  

% Verification case markers
case_colors = lines(N_cases);
for k = 1:N_cases
    St_k  = results(k, 1);
    v0_k  = results(k, 2);
    err_k = results(k, 6);

    % Circle at the (St, v0/Uref) of each case
    if k == 1
        loglog(ax1, St_k, v0_k, 'o', ...
            'Color', 'k', 'MarkerFaceColor', case_colors(k,:), 'MarkerSize', 10, ...
            'DisplayName', 'Verification cases');
    else
        loglog(ax1, St_k, v0_k, 'o', ...
            'Color', 'k', 'MarkerFaceColor', case_colors(k,:), 'MarkerSize', 10, ...
            'HandleVisibility', 'off');
    end

    % Label: case number + error
    text(ax1, St_k * 1.15, v0_k, ...
        sprintf('  %d  (%.3f%%)', k, err_k), ...
        'FontSize', 9.5, 'Color', 'k', 'FontName', 'Arial', ...
        'VerticalAlignment', 'middle');
end

%% Axes and labels
xlabel(ax1, 'St = \tau_p U_{ref} / H', 'FontSize', 13);
ylabel(ax1, 'v_0 / U_{ref}',           'FontSize', 13);
title(ax1, {'Gravity-induced rigid translation of transport regime boundaries', ...
    sprintf('Shift = \\Gamma \\cdot St,   \\Gamma = gH/U_{ref}^2 = %.4f', Gamma_outer)}, ...
    'FontSize', 12, 'FontWeight', 'normal');
grid(ax1, 'on'); box(ax1, 'on');
ax1.XMinorGrid = 'on'; ax1.YMinorGrid = 'on';
xlim(ax1, [St_plot(1), St_plot(end)]);
ylim(ax1, [v0_lo, v0_hi]);
ax1.FontSize = 12; ax1.FontName = 'Arial';

% Regime zone labels 
text(ax1, 8, 0.14, 'Near-wall confinement', ...
    'FontSize', 10, 'Color', [0.10 0.14 0.34]*0.6, 'FontName', 'Arial', ...
    'HorizontalAlignment', 'right', 'FontAngle', 'italic');
text(ax1, 8, 0.65, 'Partial penetration', ...
    'FontSize', 10, 'Color', [0.15 0.40 0.65], 'FontName', 'Arial', ...
    'HorizontalAlignment', 'right', 'FontAngle', 'italic');
text(ax1, 8, 7.0,  'Midplane crossing', ...
    'FontSize', 10, 'Color', [0.65 0.50 0.05], 'FontName', 'Arial', ...
    'HorizontalAlignment', 'right', 'FontAngle', 'italic');

% Legend
leg1 = legend(ax1, 'Location', 'southwest', 'FontSize', 10);
leg1.Box = 'on';

saveas(fh1, 'fig_gravity_translation.png');
fprintf('Saved → fig_gravity_translation.png\n');

%% FIGURE 2: Trajectory pairs (gravity off vs on)
case_colors = [0.20 0.45 0.75;   % blue
               0.85 0.20 0.10;   % red
               0.15 0.60 0.30;   % green
               0.60 0.20 0.70;   % purple
               0.90 0.55 0.00];  % orange

fh2 = figure('Name', 'Fig — Gravity trajectory pairs', ...
             'Position', [0 0 1100 680], 'Color', 'w');

for k = 1:N_cases
    tau_p = tau_p_cases(k);
    v0    = v0_cases(k);
    St_k  = tau_p * Uavg / H;
    v0_k  = v0 / Uavg;
    col   = case_colors(k,:);

    % Gravity OFF trajectory
    [~, s_off] = single_particle_run_history_gravity( ...
        tau_p, v0, 0, H, L_verif, Uavg, y0_init, ...
        p_stick_verif, N_bounce_verif, target, ...
        tmax_verif, RelTol, AbsTol, MaxStep, false, g_phys, true);

    % Gravity ON trajectory
    [~, s_on] = single_particle_run_history_gravity( ...
        tau_p, v0, 0, H, L_verif, Uavg, y0_init, ...
        p_stick_verif, N_bounce_verif, target, ...
        tmax_verif, RelTol, AbsTol, MaxStep, true, g_phys, true);

    % Normalise to outer units
    x_off = s_off(:,1) ./ H;
    y_off = 2 .* (s_off(:,2) ./ H) - 1;
    x_on  = s_on(:,1)  ./ H;
    y_on  = 2 .* (s_on(:,2)  ./ H) - 1;

    % y_max locations
    [ymax_off_val, imax_off] = max(y_off);
    [ymax_on_val,  imax_on]  = max(y_on);

    % Subplot
    ax = subplot(2, 3, k, 'Parent', fh2);
    hold(ax, 'on');

    % Channel walls and midplane
    yline(ax, -1, 'k-',  'LineWidth', 1.0, 'HandleVisibility', 'off');
    yline(ax,  0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off', 'Alpha', 0.5);

    % Gravity-off trajectory and peak
    plot(ax, x_off, y_off, '-', 'Color', col, 'LineWidth', 2.0, ...
         'DisplayName', 'No gravity');
    plot(ax, x_off(imax_off), ymax_off_val, 'o', ...
         'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 8, ...
         'HandleVisibility', 'off');

    % Gravity-on trajectory and peak
    plot(ax, x_on, y_on, '--', 'Color', col*0.65, 'LineWidth', 2.0, ...
         'DisplayName', 'With gravity');
    plot(ax, x_on(imax_on), ymax_on_val, 's', ...
         'Color', col*0.65, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col*0.65, ...
         'MarkerSize', 8, 'HandleVisibility', 'off');

    % Horizontal dotted lines at the two y_max values
    yline(ax, ymax_off_val, ':', 'Color', col,      'LineWidth', 1.2, 'Alpha', 0.7, ...
          'HandleVisibility', 'off');
    yline(ax, ymax_on_val,  ':', 'Color', col*0.65, 'LineWidth', 1.2, 'Alpha', 0.7, ...
          'HandleVisibility', 'off');

    % Gravity shift arrow: from ymax_on up to ymax_off at x midpoint
    x_mid_arr = mean(xlim(ax));
    if abs(ymax_off_val - ymax_on_val) > 0.005
        annotation_x = 0.5;  
        text(ax, x_off(imax_off)*0.6, (ymax_off_val + ymax_on_val)/2, ...
            sprintf('\\Delta y_{max}\n= %.4f', ymax_off_val - ymax_on_val), ...
            'FontSize', 8.5, 'Color', [0.4 0.4 0.4], 'HorizontalAlignment', 'center');
    end

    % Zoom y-axis 
    y_range_lo = min(-1.02, ymax_on_val  - 0.05);
    y_range_hi = max(ymax_off_val + 0.08, -0.80);
    ylim(ax, [y_range_lo, y_range_hi]);

    % Axes and title
    xlabel(ax, 'x / H', 'FontSize', 10);
    if mod(k,3) == 1
        ylabel(ax, 'y  (outer units)', 'FontSize', 10);
    end
    title(ax, sprintf('Case %d:  St = %.3f,  v_0/U = %.2f\nAnalytical shift = %.5f,  Error = %.4f%%', ...
        k, St_k, v0_k, results(k,5) - results(k,3), results(k,6)), ...
        'FontSize', 9.5, 'FontWeight', 'normal');

    grid(ax, 'on'); box(ax, 'on');
    ax.FontSize = 10; ax.FontName = 'Arial';

    if k == 1
        legend(ax, 'Location', 'southeast', 'FontSize', 9);
    end
end

% Global title
sgtitle(fh2, sprintf(['Gravity verification trajectories — no gravity (solid) vs gravity (dashed)\n' ...
    '\\Gamma = gH/U_{ref}^2 = %.4f.   Filled circle = y_{max} no grav;   Open square = y_{max} with grav.'], ...
    Gamma_outer), 'FontSize', 11, 'FontName', 'Arial');

ax6 = subplot(2,3,6,'Parent',fh2);
set(ax6,'Visible','off');

axes_pos = get(ax6, 'Position');
str_lines = {'Summary of relative errors:'; ''};
for k = 1:N_cases
    str_lines{end+1} = sprintf('Case %d  St=%.3f  v_0/U=%.2f  err=%.4f%%', ...
        k, results(k,1), results(k,2), results(k,6)); %#ok<AGROW>
end
str_lines{end+1} = '';
str_lines{end+1} = sprintf('Max error: %.4f%%', max_error);
if max_error < 2.0
    str_lines{end+1} = 'R6: SATISFIED (< 2%)';
else
    str_lines{end+1} = 'R6: NOT SATISFIED (> 2%)';
end
annotation(fh2, 'textbox', axes_pos, ...
    'String', str_lines, 'EdgeColor', [0.8 0.8 0.8], ...
    'BackgroundColor', [0.97 0.97 0.97], ...
    'FitBoxToText', 'off', 'FontSize', 9.5, 'FontName', 'Arial', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

saveas(fh2, 'fig_gravity_trajectories.png');
fprintf('Saved → fig_gravity_trajectories.png\n\n');

fprintf('=============================================================\n');
fprintf('  All figures and results saved successfully.\n');
fprintf('=============================================================\n');

end  % end main function


%% ── Single particle run with optional gravity (scalar outputs) ───────────────

function [outcome, x_final, y_final, t_total, crossed_target, ...
          bounce_at_stick, y_max] = single_particle_run_gravity( ...
          tau_p, v0_perp, x0, ...
          H, L, Uavg, y0, ...
          p_stick, N_bounce_max, target, ...
          tmax, RelTol, AbsTol, MaxStep, ...
          useGravity, g_phys, useUniformFlow)

if nargin < 18, useUniformFlow = false; end
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
    rhs     = @(t,s) particle_rhs_gravity(s, tau_p, H, Uavg, useGravity, g_phys, useUniformFlow);
    eventsf = @(t,s) particle_events_grav(s, H, L, target);
    opts    = odeset('RelTol', RelTol, 'AbsTol', AbsTol, ...
                     'MaxStep', MaxStep, 'Events', eventsf);

    [t_sol, s_sol, ~, se, ie] = ode45(rhs, [0, tmax - t_total], s0, opts);

    t_total = t_total + t_sol(end);
    y_max   = max(y_max, max(s_sol(:,2)));

    if isempty(ie)
        outcome = 0;
        x_final = s_sol(end,1); y_final = s_sol(end,2);
        break;
    end

    if any(ie == 4), crossed_target = true; end

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

    % Wall contact: stochastic stick/bounce — same as laminar model
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

    % Elastic reflection — identical offset and clamp as laminar model
    xR = sv(1); yR = sv(2); vxR = sv(3); vyR = sv(4);
    if ev == 1
        vyR = abs(vyR);   yR = max(yR, 0) + 1e-12;   % bottom wall bounce
    else
        vyR = -abs(vyR);  yR = min(yR, H) - 1e-12;   % top wall bounce
    end
    s0 = [xR; yR; vxR; vyR];
end

end


%% Single particle run with full trajectory history

function [t_hist, s_hist] = single_particle_run_history_gravity( ...
          tau_p, v0_perp, x0, ...
          H, L, Uavg, y0, ...
          p_stick, N_bounce_max, target, ...
          tmax, RelTol, AbsTol, MaxStep, ...
          useGravity, g_phys, useUniformFlow)

if nargin < 18, useUniformFlow = false; end

s0      = [x0; y0; 0.0; v0_perp];
bounced = 0;
t_total = 0.0;
t_segs  = {};
s_segs  = {};

while true
    rhs     = @(t,s) particle_rhs_gravity(s, tau_p, H, Uavg, useGravity, g_phys, useUniformFlow);
    eventsf = @(t,s) particle_events_grav(s, H, L, target);
    opts    = odeset('RelTol', RelTol, 'AbsTol', AbsTol, ...
                     'MaxStep', MaxStep, 'Events', eventsf);

    [t_sol, s_sol, ~, se, ie] = ode45(rhs, [0, tmax - t_total], s0, opts);

    if isempty(t_segs)
        t_segs{end+1} = t_sol + t_total;        %#ok<AGROW>
        s_segs{end+1} = s_sol;                   %#ok<AGROW>
    else
        t_segs{end+1} = t_sol(2:end) + t_total; %#ok<AGROW>
        s_segs{end+1} = s_sol(2:end,:);          %#ok<AGROW>
    end
    t_total = t_total + t_sol(end);

    if isempty(ie), break; end

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
        vyR = abs(vyR);   yR = max(yR, 0) + 1e-12;
    else
        vyR = -abs(vyR);  yR = min(yR, H) - 1e-12;
    end
    s0 = [xR; yR; vxR; vyR];
end

t_hist = vertcat(t_segs{:});
s_hist = vertcat(s_segs{:});

end


%% ODE RHS: Stokes drag + optional gravity + optional uniform streamwise flow 

function ds = particle_rhs_gravity(s, tau_p, H, Uavg, useGravity, g_phys, useUniformFlow)

if nargin < 7, useUniformFlow = false; end

y  = s(2);
vx = s(3);
vy = s(4);

if useUniformFlow
    u_f = Uavg;          % constant streamwise flow 
else
    yclamp = min(max(y, 0), H);
    u_f    = 6 * Uavg * (yclamp/H) * (1 - yclamp/H);
end
v_f = 0;

ax = (u_f - vx) / tau_p;
ay = (v_f - vy) / tau_p;

if useGravity
    ay = ay - g_phys;
end

ds = [vx; vy; ax; ay];

end


%% Events: bottom wall (ev=1), top wall (ev=2), outlet (ev=3), target (ev=4)
function [value, isterminal, direction] = particle_events_grav(s, H, L, target)

x = s(1); y = s(2);

v_bottom = y;
v_top    = H - y;
v_out    = L - x;

if target.enabled
    dx    = x - target.x; dy = y - target.y;
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