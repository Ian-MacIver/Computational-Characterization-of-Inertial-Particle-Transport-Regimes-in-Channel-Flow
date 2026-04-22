function compute_ensemble_statistics()
% Input files:  turbulent_particle_results_#.mat
% Output file:  turbulent_ensemble_results.mat

clearvars;
clc;

%% Settings
N_SNAPSHOTS    = 11;
RESULT_PREFIX  = 'turbulent_particle_final_results_'; % Ensure # match
ENSEMBLE_FILE  = 'turbulent_ensemble_results.mat';

% Regime boundary thresholds
THRESH_PARTIAL = -0.5;   % y_max threshold for partial penetration boundary
THRESH_MID     =  0.0;   % y_max threshold for midplane crossing boundary

% Minimum valid points required for a power-law fit at each boundary
MIN_FIT_POINTS = 4;

fig_dir = 'figures_ensemble';
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

fprintf('=== Ensemble aggregation: %d snapshots ===\n', N_SNAPSHOTS);

%% Load and validate all snapshot files
ref_file = sprintf('%s%d.mat', RESULT_PREFIX, 1);
if ~isfile(ref_file)
    error('Cannot find %s. Run all eleven sweeps first.', ref_file);
end

R = load(ref_file, 'tau_p_list', 'v0_list', 'St_list', 'y_max_map', 'params');
tau_p_list = R.tau_p_list;
v0_list    = R.v0_list;
St_list    = R.St_list;
ref_params = R.params;

Ntau = numel(tau_p_list);
Nv0  = numel(v0_list);

% Stack array: [Ntau x Nv0 x N_SNAPSHOTS]
y_max_stack       = nan(Ntau, Nv0, N_SNAPSHOTS);
outcome_stack     = nan(Ntau, Nv0, N_SNAPSHOTS);
t_res_stack       = nan(Ntau, Nv0, N_SNAPSHOTS);
bounce_stack      = nan(Ntau, Nv0, N_SNAPSHOTS);
snapshot_runtimes = nan(N_SNAPSHOTS, 1);

for k = 1:N_SNAPSHOTS
    fname = sprintf('%s%d.mat', RESULT_PREFIX, k);

    if ~isfile(fname)
        error('Missing result file: %s\nRun snapshot %d before aggregating.', fname, k);
    end

    S = load(fname, 'tau_p_list', 'v0_list', 'St_list', ...
                    'y_max_map', 'outcome_map', 't_res_map', ...
                    'bounce_at_stick_map', 'params', 'diagnostics');

    % Validate parameter consistency across snapshots
    if ~isequal(S.tau_p_list, tau_p_list) || ~isequal(S.v0_list, v0_list)
        error('Parameter vectors in %s do not match snapshot 1. Re-run sweep %d.', fname, k);
    end
    if ~isequal(S.params.solver, ref_params.solver) || ...
       ~isequal(S.params.wall,   ref_params.wall)   || ...
       ~isequal(S.params.dom,    ref_params.dom)
        warning('Solver/wall/domain settings differ in snapshot %d — results may not be comparable.', k);
    end

    % Squeeze out the singleton x0 dimension (Ntau x Nv0 x 1 -> Ntau x Nv0)
    y_max_stack(:,:,k)   = squeeze(S.y_max_map(:,:,1));
    outcome_stack(:,:,k) = squeeze(double(S.outcome_map(:,:,1)));
    t_res_stack(:,:,k)   = squeeze(S.t_res_map(:,:,1));
    bounce_stack(:,:,k)  = squeeze(S.bounce_at_stick_map(:,:,1));
    snapshot_runtimes(k) = S.diagnostics.runtime_sec;

    fprintf('  Loaded snapshot %2d  (runtime: %.1f min)\n', k, snapshot_runtimes(k)/60);
end

fprintf('All %d snapshots loaded successfully.\n\n', N_SNAPSHOTS);

%% Ensemble statistics
% Pointwise mean and standard deviation across snapshots
ensemble_mean_ymax = mean(y_max_stack,    3, 'omitnan');
ensemble_std_ymax  = std(y_max_stack,  0, 3, 'omitnan');
ensemble_mean_outcome = mean(outcome_stack, 3, 'omitnan');
ensemble_mean_tres    = mean(t_res_stack,   3, 'omitnan');

fprintf('Ensemble mean y_max range:  [%.3f, %.3f]\n', ...
    min(ensemble_mean_ymax(:)), max(ensemble_mean_ymax(:)));
fprintf('Ensemble std  y_max range:  [%.4f, %.4f]\n', ...
    min(ensemble_std_ymax(:)),  max(ensemble_std_ymax(:)));

%% Per-snapshot boundary fitting
% For each snapshot, extract v0_crit(St) at both thresholds and fit
% v0_crit = C * St^alpha by linear regression in log-log space.

ref_Uref = ref_params.ref.Uref;
v0_norm  = v0_list ./ ref_Uref;

alpha_partial_vec = nan(N_SNAPSHOTS, 1);
alpha_mid_vec     = nan(N_SNAPSHOTS, 1);
C_partial_vec     = nan(N_SNAPSHOTS, 1);
C_mid_vec         = nan(N_SNAPSHOTS, 1);

for k = 1:N_SNAPSHOTS
    ymap = y_max_stack(:,:,k);

    [alpha_partial_vec(k), C_partial_vec(k)] = ...
        fit_power_law_boundary(St_list, v0_norm, ymap, THRESH_PARTIAL, MIN_FIT_POINTS);
    [alpha_mid_vec(k), C_mid_vec(k)] = ...
        fit_power_law_boundary(St_list, v0_norm, ymap, THRESH_MID, MIN_FIT_POINTS);
end

fprintf('\nPer-snapshot boundary fits:\n');
fprintf('  %-10s  %-12s  %-12s  %-12s  %-12s\n', ...
    'Snapshot', 'alpha_part', 'C_part', 'alpha_mid', 'C_mid');
for k = 1:N_SNAPSHOTS
    fprintf('  %-10d  %-12.4f  %-12.4f  %-12.4f  %-12.4f\n', ...
        k, alpha_partial_vec(k), C_partial_vec(k), ...
        alpha_mid_vec(k), C_mid_vec(k));
end

%% Ensemble boundary fit (from mean map)
[alpha_partial_ensemble, C_partial_ensemble] = ...
    fit_power_law_boundary(St_list, v0_norm, ensemble_mean_ymax, THRESH_PARTIAL, MIN_FIT_POINTS);
[alpha_mid_ensemble, C_mid_ensemble] = ...
    fit_power_law_boundary(St_list, v0_norm, ensemble_mean_ymax, THRESH_MID, MIN_FIT_POINTS);

fprintf('\nEnsemble mean map boundary fits:\n');
fprintf('  Partial penetration:  alpha = %.4f,  C = %.4f\n', alpha_partial_ensemble, C_partial_ensemble);
fprintf('  Midplane crossing:    alpha = %.4f,  C = %.4f\n', alpha_mid_ensemble,     C_mid_ensemble);

%% Convergence assessment
% Running mean of alpha_mid as snapshots are added sequentially.
% Used to confirm that eleven snapshots are statistically sufficient.
running_mean_alpha_partial = nan(N_SNAPSHOTS, 1);
running_mean_alpha_mid     = nan(N_SNAPSHOTS, 1);

for n = 1:N_SNAPSHOTS
    running_mean_alpha_partial(n) = mean(alpha_partial_vec(1:n), 'omitnan');
    running_mean_alpha_mid(n)     = mean(alpha_mid_vec(1:n),     'omitnan');
end

% Find the first snapshot at which running mean is within 5% of final value
final_alpha_mid = running_mean_alpha_mid(end);
converged_at = find(abs(running_mean_alpha_mid - final_alpha_mid) ...
                    ./ abs(final_alpha_mid) < 0.05, 1, 'first');

fprintf('\nConvergence assessment (midplane boundary exponent alpha_mid):\n');
fprintf('  Final ensemble value: %.4f\n', final_alpha_mid);
if ~isempty(converged_at)
    fprintf('  Running mean converges within 5%% of final value at snapshot %d\n', converged_at);
else
    fprintf('  WARNING: Running mean did not converge within 5%% — consider additional snapshots.\n');
end

%% Save ensemble results
fprintf('\nSaving ensemble results to %s ...\n', ENSEMBLE_FILE);

save(ENSEMBLE_FILE, ...
    'tau_p_list', 'v0_list', 'St_list', 'v0_norm', ...
    'y_max_stack', 'outcome_stack', 't_res_stack', 'bounce_stack', ...
    'ensemble_mean_ymax', 'ensemble_std_ymax', ...
    'ensemble_mean_outcome', 'ensemble_mean_tres', ...
    'alpha_partial_vec', 'C_partial_vec', ...
    'alpha_mid_vec',     'C_mid_vec', ...
    'alpha_partial_ensemble', 'C_partial_ensemble', ...
    'alpha_mid_ensemble',     'C_mid_ensemble', ...
    'running_mean_alpha_partial', 'running_mean_alpha_mid', ...
    'converged_at', 'snapshot_runtimes', ...
    'N_SNAPSHOTS', 'THRESH_PARTIAL', 'THRESH_MID', ...
    'ref_params', '-v7.3');

fprintf('Done. Ensemble results saved.\n');
fprintf('\n=== Summary for thesis ===\n');
fprintf('  Ensemble alpha_mid:     %.4f  (laminar prediction: -1.000)\n', alpha_mid_ensemble);
fprintf('  Ensemble alpha_partial: %.4f\n', alpha_partial_ensemble);
fprintf('  Convergence snapshot:   %d of %d\n', converged_at, N_SNAPSHOTS);
fprintf('  Total sweep time:       %.1f h\n', sum(snapshot_runtimes)/3600);
fprintf('  Mean per-snapshot time: %.1f min\n', mean(snapshot_runtimes)/60);

end


%% Local helper: power-law boundary fit
function [alpha, C] = fit_power_law_boundary(St_list, v0_norm, ymax_map, threshold, min_pts)
% For each St, find the lowest v0 at which ymax exceeds threshold.
% Fit log(v0_crit) = alpha*log(St) + log(C) by linear regression.
% Returns NaN if fewer than min_pts valid points are found.

Ntau = numel(St_list);
v0_thresh = nan(Ntau, 1);

for ia = 1:Ntau
    row  = ymax_map(ia, :);
    idx  = find(row > threshold, 1, 'first');
    if ~isempty(idx)
        v0_thresh(ia) = v0_norm(idx);
    end
end

valid = ~isnan(v0_thresh) & v0_thresh > 0 & St_list(:) > 0;

if sum(valid) < min_pts
    alpha = NaN;
    C     = NaN;
    return;
end

p     = polyfit(log(St_list(valid)), log(v0_thresh(valid)), 1);
alpha = p(1);
C     = exp(p(2));

end
