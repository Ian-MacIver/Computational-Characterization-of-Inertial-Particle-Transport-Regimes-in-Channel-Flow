function run_turbulent_particle_sweep_Final()
%% Housekeeping
clearvars;
clc;

set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLineLineWidth',1.2);
set(0,'DefaultFigureColor','w');

rng('default');

run_name        = 'channel5200_frozen_snapshot_run';
results_file    = 'turbulent_particle_final_results_1.mat'; % edit # to create new file
flow_cache_file = 'channel5200_flow_cache_11.mat'; % edit # to create new file

%% Channel Geometry
geom.y_bottom = -1.0;
geom.y_top    =  1.0;
geom.H        = geom.y_top - geom.y_bottom;

%% Reference Scales
ref.Uref     = 1.0;
ref.St_label = 'St = \tau_p U_{ref}/H';

%% Turbulence Config
turb.useTurbulence    = true;
turb.authkey          = 'auth-key'; % request from JHTDB
turb.dataset          = 'channel5200';
turb.var              = 'velocity';
turb.temporal_method  = 'none';    % frozen snapshot
turb.spatial_method   = 'lag6';
turb.spatial_operator = 'field';
turb.time0            = 1;         % valid integer snapshot in [1,11]
turb.Lx               = 8*pi;
turb.Lz               = 3*pi;

if turb.time0 < 1 || turb.time0 > 11 || turb.time0 ~= round(turb.time0)
    error('For channel5200, turb.time0 must be an integer from 1 to 11.');
end

%% Force Model (gravity)
forces.useGravity = false;
forces.g          = [0; -9.81; 0];

%% Interpolation Domain
dom.x_min =  0.0;
dom.x_max =  8*pi;   % full periodic length of channel5200

dom.y_min = -1.0;
dom.y_max =  1.0;

dom.z_min = -1.0;
dom.z_max =  1.0;

dom.Nx = 96;
dom.Ny = 72;
dom.Nz = 64;

dom.x = linspace(dom.x_min, dom.x_max, dom.Nx);
dom.y = linspace(dom.y_min, dom.y_max, dom.Ny);
dom.z = linspace(dom.z_min, dom.z_max, dom.Nz);

%% Sweep Definitions

% Refined final sweep
tau_p_list = logspace(-3, 1, 41);
v0_list    = linspace(0.1*ref.Uref, 6.0*ref.Uref, 41);
x0_list    = 0.0;
z0         = 0.0;

St_list = tau_p_list .* (ref.Uref ./ geom.H);

Ntau = numel(tau_p_list);
Nv0  = numel(v0_list);
Nx0  = numel(x0_list);

%% Particle-wall Model
wall.p_stick_wall = 0.20;
wall.N_bounce_max = 50;

%% Target Options
target.enabled = false;
target.x       = 1.0;
target.y       = 0.0;
target.z       = 0.0;
target.eps     = 0.05;

%% Solver Settings


% Thesis settings
solver.name        = 'adaptive';   % ode45 for tau_p >= tau_p_stiff, ode15s otherwise
solver.tmax    = 50.0;
solver.RelTol  = 1e-6;
solver.AbsTol  = 1e-9;
solver.MaxStep = 0.02;
solver.tau_p_stiff = 0.1;          % stiffness threshold in outer units


%% Flow Field
if isfile(flow_cache_file)
    S = load(flow_cache_file, 'flow_full', 'dom_saved', 'turb_saved');
    if isequal(S.dom_saved, dom) && isequaln(S.turb_saved, turb)
        flow_full = S.flow_full;
        fprintf('Loaded cached flow from %s\n', flow_cache_file);
    else
        fprintf('Cached flow does not match current settings. Recomputing...\n');
        flow_full  = precompute_velocity_domain(turb, dom);
        dom_saved  = dom;
        turb_saved = turb;
        save(flow_cache_file, 'flow_full', 'dom_saved', 'turb_saved', '-v7.3');
    end
else
    fprintf('Precomputing frozen velocity snapshot on %d x %d x %d domain...\n', ...
        dom.Nx, dom.Ny, dom.Nz);
    flow_full  = precompute_velocity_domain(turb, dom);
    dom_saved  = dom;
    turb_saved = turb;
    save(flow_cache_file, 'flow_full', 'dom_saved', 'turb_saved', '-v7.3');
end
fprintf('Velocity field ready.\n');

% Slim flow struct for parfor broadcast (no raw U/V/W arrays).
flow_slim.dom = flow_full.dom;
flow_slim.Fu  = flow_full.Fu;
flow_slim.Fv  = flow_full.Fv;
flow_slim.Fw  = flow_full.Fw;

%% Flatten Sweep
[TauGrid, V0Grid, X0Grid] = ndgrid(tau_p_list, v0_list, x0_list);
tau_cases = TauGrid(:);
v0_cases  = V0Grid(:);
x0_cases  = X0Grid(:);

Ncases = numel(tau_cases);

outcome_vec         = zeros(Ncases, 1);
final_x_vec         = nan(Ncases, 1);
final_y_vec         = nan(Ncases, 1);
final_z_vec         = nan(Ncases, 1);
crossed_target_vec  = false(Ncases, 1);
bounce_at_stick_vec = nan(Ncases, 1);
t_res_vec           = nan(Ncases, 1);
y_max_vec           = nan(Ncases, 1);

%% Parallel Pool
pool = gcp('nocreate');
if isempty(pool)
    parpool;
end

fprintf('Starting parallel sweep with %d cases...\n', Ncases);
tic;

parfor idx = 1:Ncases
    tau_p   = tau_cases(idx);
    v0_perp = v0_cases(idx);
    x0      = x0_cases(idx);

    [outcome, x_final, y_final, z_final, t_total, crossed_target, ...
        bounce_at_stick, y_max] = single_particle_run_3d( ...
        tau_p, v0_perp, x0, z0, ...
        geom.y_bottom, geom.y_top, ...
        wall.p_stick_wall, wall.N_bounce_max, target, ...
        flow_slim, forces, solver);

    outcome_vec(idx)         = outcome;
    final_x_vec(idx)         = x_final;
    final_y_vec(idx)         = y_final;
    final_z_vec(idx)         = z_final;
    crossed_target_vec(idx)  = crossed_target;
    bounce_at_stick_vec(idx) = bounce_at_stick;
    t_res_vec(idx)           = t_total;
    y_max_vec(idx)           = y_max;
end

runtime_sec = toc;
fprintf('Parallel sweep complete in %.1f s.\n', runtime_sec);

%% Reshape
outcome_map         = reshape(outcome_vec,         [Ntau, Nv0, Nx0]);
final_x_map         = reshape(final_x_vec,         [Ntau, Nv0, Nx0]);
final_y_map         = reshape(final_y_vec,         [Ntau, Nv0, Nx0]);
final_z_map         = reshape(final_z_vec,         [Ntau, Nv0, Nx0]);
crossed_target_map  = reshape(crossed_target_vec,  [Ntau, Nv0, Nx0]);
bounce_at_stick_map = reshape(bounce_at_stick_vec, [Ntau, Nv0, Nx0]);
t_res_map           = reshape(t_res_vec,           [Ntau, Nv0, Nx0]);
y_max_map           = reshape(y_max_vec,           [Ntau, Nv0, Nx0]);

%% Representative Trajectories
% St_list(6)  ~ 0.003  (tracer)
% St_list(13) ~ 0.1    (near transition - midpoint of 25-pt logspace(-3,1))
% St_list(20) ~ 3.0    (ballistic)
% v0_list(ceil(Nv0*0.4)) gives a moderate injection speed

traj_cases = struct([]);

traj_cases(1).label = 'low_St_near_wall';
traj_cases(1).tau_p = tau_p_list(6);
traj_cases(1).v0    = v0_list(ceil(Nv0 * 0.6));   % moderate v0

traj_cases(2).label = 'mid_St_transitional';
traj_cases(2).tau_p = tau_p_list(13);              % near St_crit
traj_cases(2).v0    = v0_list(ceil(Nv0 * 0.6));

traj_cases(3).label = 'high_St_ballistic';
traj_cases(3).tau_p = tau_p_list(20);              % well into ballistic regime
traj_cases(3).v0    = v0_list(ceil(Nv0 * 0.4));   % moderate v0 -> visible bouncing

trajectories = struct([]);
for k = 1:numel(traj_cases)
    [t_hist, s_hist, info] = single_particle_run_3d_history( ...
        traj_cases(k).tau_p, traj_cases(k).v0, 0.0, z0, ...
        geom.y_bottom, geom.y_top, ...
        wall.p_stick_wall, wall.N_bounce_max, target, ...
        flow_slim, forces, solver);

    trajectories(k).label = traj_cases(k).label;
    trajectories(k).tau_p = traj_cases(k).tau_p;
    trajectories(k).v0    = traj_cases(k).v0;
    trajectories(k).t     = t_hist;
    trajectories(k).s     = s_hist;
    trajectories(k).info  = info;
end

%% Diagnostics
diagnostics.frac_hit_tmax = mean(t_res_map(:) >= solver.tmax);
diagnostics.near_x_edge   = mean(final_x_map(:) <= dom.x_min + 1e-3 | ...
                                  final_x_map(:) >= dom.x_max - 1e-3, 'omitnan');
diagnostics.near_z_edge   = mean(final_z_map(:) <= dom.z_min + 1e-3 | ...
                                  final_z_map(:) >= dom.z_max - 1e-3, 'omitnan');
diagnostics.frac_midplane = mean(y_max_map(:) > 0.0, 'omitnan');
diagnostics.runtime_sec   = runtime_sec;
diagnostics.num_cases     = Ncases;

fprintf('Fraction hitting time cap:    %.3f\n', diagnostics.frac_hit_tmax);
fprintf('Fraction ending near x-edge:  %.3f\n', diagnostics.near_x_edge);
fprintf('Fraction ending near z-edge:  %.3f\n', diagnostics.near_z_edge);
fprintf('Fraction with y_max > 0:      %.3f\n', diagnostics.frac_midplane);

%% Save
params.run_name = run_name;
params.geom     = geom;
params.ref      = ref;
params.turb     = turb;
params.forces   = forces;
params.dom      = dom;
params.wall     = wall;
params.target   = target;
params.solver   = solver;

save(results_file, ...
    'run_name', ...
    'tau_p_list','v0_list','x0_list','z0','St_list', ...
    'outcome_map','final_x_map','final_y_map','final_z_map', ...
    'crossed_target_map','bounce_at_stick_map','t_res_map','y_max_map', ...
    'trajectories','diagnostics','params', ...
    '-v7.3');

fprintf('Saved results to %s\n', results_file);

end


%% Precompute frozen 3D velocity domain
function flow_full = precompute_velocity_domain(turb, dom)

Nx = dom.Nx;
Ny = dom.Ny;
Nz = dom.Nz;
N  = Nx * Ny * Nz;

% Vectorised point array construction via ndgrid (replaces triple loop).
[Xd, Yd, Zd] = ndgrid(dom.x, dom.y, dom.z);
points = [Xd(:), Yd(:), Zd(:)];   % N x 3

result = JHUTDB_getData( ...
    turb.authkey, turb.dataset, turb.var, turb.time0, ...
    turb.temporal_method, turb.spatial_method, turb.spatial_operator, ...
    points);

if ~isnumeric(result) || size(result,2) ~= 3 || size(result,1) ~= N
    error('Velocity query did not return an N x 3 numeric array.');
end

% ndgrid ordering -> reshape directly to [Nx, Ny, Nz].
U = reshape(result(:,1), [Nx, Ny, Nz]);
V = reshape(result(:,2), [Nx, Ny, Nz]);
W = reshape(result(:,3), [Nx, Ny, Nz]);

flow_full.dom = dom;
flow_full.U   = U;   % kept for caching / post-processing only
flow_full.V   = V;
flow_full.W   = W;

flow_full.Fu = griddedInterpolant(Xd, Yd, Zd, U, 'linear', 'nearest');
flow_full.Fv = griddedInterpolant(Xd, Yd, Zd, V, 'linear', 'nearest');
flow_full.Fw = griddedInterpolant(Xd, Yd, Zd, W, 'linear', 'nearest');

end


%% One Particle Trajectory (sweep)
function [outcome, x_final, y_final, z_final, t_total, crossed_target, ...
    bounce_at_stick, y_max] = single_particle_run_3d( ...
    tau_p, v0_perp, x0, z0, ...
    y_bottom, y_top, ...
    p_stick_wall, N_bounce_max, target, ...
    flow, forces, solver)

outcome         = 0;
x_final         = NaN;
y_final         = NaN;
z_final         = NaN;
t_total         = 0.0;
crossed_target  = false;
bounce_at_stick = NaN;

eps_wall = 1e-6;

x  = x0;
y  = y_bottom + eps_wall;
z  = z0;
vx = 0.0;
vy = v0_perp;
vz = 0.0;

s0    = [x; y; z; vx; vy; vz];
y_max = y;
bounced = 0;

% Build ODE options and function handles once outside the bounce loop.
rhs     = @(t,s) particle_rhs_3d_interp(t, s, tau_p, flow, forces);
eventsf = @(t,s) particle_events_3d_interp(t, s, y_bottom, y_top, target, flow.dom);
opts    = odeset('RelTol', solver.RelTol, 'AbsTol', solver.AbsTol, ...
                 'MaxStep', solver.MaxStep, 'Events', eventsf);

% Adaptive solver: stiff (ode15s) for small tau_p, non-stiff (ode45) otherwise.
if tau_p < solver.tau_p_stiff
    ode_solver = @ode15s;
else
    ode_solver = @ode45;
end

while true
    [t_sol, s_sol, ~, se, ie] = ode_solver(rhs, [0, solver.tmax - t_total], s0, opts);

    t_total = t_total + t_sol(end);
    y_max   = max(y_max, max(s_sol(:,2)));

    if isempty(ie)
        outcome = 0;
        x_final = s_sol(end,1);
        y_final = s_sol(end,2);
        z_final = s_sol(end,3);
        break;
    end

    if any(ie == 3)
        crossed_target = true;
    end

    term_mask = (ie == 1) | (ie == 2) | (ie == 4);

    if ~any(term_mask)
        if t_total >= solver.tmax
            outcome = 0;
            x_final = s_sol(end,1);
            y_final = s_sol(end,2);
            z_final = s_sol(end,3);
            break;
        else
            s0 = s_sol(end,:).';
            continue;
        end
    end

    idx_last = find(term_mask, 1, 'last');
    ev = ie(idx_last);
    sv = se(idx_last,:);

    if ev == 4
        outcome = 0;
        x_final = sv(1);
        y_final = sv(2);
        z_final = sv(3);
        break;
    end

    if rand() <= p_stick_wall
        outcome         = 1;
        x_final         = sv(1);
        y_final         = sv(2);
        z_final         = sv(3);
        bounce_at_stick = bounced;
        break;
    end

    bounced = bounced + 1;

    if bounced >= N_bounce_max
        outcome         = 1;
        x_final         = sv(1);
        y_final         = sv(2);
        z_final         = sv(3);
        bounce_at_stick = bounced;
        break;
    end

    xR  = sv(1);
    yR  = sv(2);
    zR  = sv(3);
    vxR = sv(4);
    vyR = sv(5);
    vzR = sv(6);

    if ev == 1
        vyR = abs(vyR);
        yR  = y_bottom + eps_wall;
    else
        vyR = -abs(vyR);
        yR  = y_top - eps_wall;
    end

    yR = max(yR, y_bottom + eps_wall);
    yR = min(yR, y_top    - eps_wall);

    s0 = [xR; yR; zR; vxR; vyR; vzR];
end

end


%% One particle trajectory with history
function [t_hist, s_hist, info] = single_particle_run_3d_history( ...
    tau_p, v0_perp, x0, z0, ...
    y_bottom, y_top, ...
    p_stick_wall, N_bounce_max, target, ...
    flow, forces, solver)

eps_wall = 1e-6;

x  = x0;
y  = y_bottom + eps_wall;
z  = z0;
vx = 0.0;
vy = v0_perp;
vz = 0.0;

s0      = [x; y; z; vx; vy; vz];
bounced = 0;
t_total = 0.0;

info   = struct('outcome', NaN, 'bounces', 0);
t_segs = {};
s_segs = {};

rhs     = @(t,s) particle_rhs_3d_interp(t, s, tau_p, flow, forces);
eventsf = @(t,s) particle_events_3d_interp(t, s, y_bottom, y_top, target, flow.dom);
opts    = odeset('RelTol', solver.RelTol, 'AbsTol', solver.AbsTol, ...
                 'MaxStep', solver.MaxStep, 'Events', eventsf);

if tau_p < solver.tau_p_stiff
    ode_solver = @ode15s;
else
    ode_solver = @ode45;
end

while true
    [t_sol, s_sol, ~, se, ie] = ode_solver(rhs, [0, solver.tmax - t_total], s0, opts);

    if isempty(t_segs)
        t_segs{end+1} = t_sol + t_total;         %#ok<AGROW>
        s_segs{end+1} = s_sol;                    %#ok<AGROW>
    else
        t_segs{end+1} = t_sol(2:end) + t_total;  %#ok<AGROW>
        s_segs{end+1} = s_sol(2:end,:);           %#ok<AGROW>
    end

    t_total = t_total + t_sol(end);

    if isempty(ie)
        info.outcome = 0;
        info.bounces = bounced;
        break;
    end

    term_mask = (ie == 1) | (ie == 2) | (ie == 4);
    if ~any(term_mask)
        if t_total >= solver.tmax
            info.outcome = 0;
            info.bounces = bounced;
            break;
        else
            s0 = s_sol(end,:).';
            continue;
        end
    end

    idx_last = find(term_mask, 1, 'last');
    ev = ie(idx_last);
    sv = se(idx_last,:);

    if ev == 4
        info.outcome = 0;
        info.bounces = bounced;
        break;
    end

    if rand() <= p_stick_wall
        info.outcome = 1;
        info.bounces = bounced;
        break;
    end

    bounced = bounced + 1;
    if bounced >= N_bounce_max
        info.outcome = 1;
        info.bounces = bounced;
        break;
    end

    xR  = sv(1);
    yR  = sv(2);
    zR  = sv(3);
    vxR = sv(4);
    vyR = sv(5);
    vzR = sv(6);

    if ev == 1
        vyR = abs(vyR);
        yR  = y_bottom + eps_wall;
    else
        vyR = -abs(vyR);
        yR  = y_top - eps_wall;
    end

    yR = max(yR, y_bottom + eps_wall);
    yR = min(yR, y_top    - eps_wall);

    s0 = [xR; yR; zR; vxR; vyR; vzR];
end

t_hist = vertcat(t_segs{:});
s_hist = vertcat(s_segs{:});

end


%% RHS
function ds = particle_rhs_3d_interp(~, s, tau_p, flow, forces)

x  = s(1);
y  = s(2);
z  = s(3);
vx = s(4);
vy = s(5);
vz = s(6);

xq = min(max(x, flow.dom.x_min), flow.dom.x_max);
yq = min(max(y, flow.dom.y_min), flow.dom.y_max);
zq = min(max(z, flow.dom.z_min), flow.dom.z_max);

uf = flow.Fu(xq, yq, zq);
vf = flow.Fv(xq, yq, zq);
wf = flow.Fw(xq, yq, zq);

ax = (uf - vx) / tau_p;
ay = (vf - vy) / tau_p;
az = (wf - vz) / tau_p;

if forces.useGravity
    ax = ax + forces.g(1);
    ay = ay + forces.g(2);
    az = az + forces.g(3);
end

ds = [vx; vy; vz; ax; ay; az];

end


%% Events
function [value, isterminal, direction] = particle_events_3d_interp(~, s, y_bottom, y_top, target, dom)

x = s(1);
y = s(2);
z = s(3);

v_bottom = y - y_bottom;
v_top    = y_top - y;

if target.enabled
    dx = x - target.x;
    dy = y - target.y;
    dz = z - target.z;
    v_target = sqrt(dx*dx + dy*dy + dz*dz) - target.eps;
else
    v_target = 1.0;
end

d_left  = x - dom.x_min;
d_right = dom.x_max - x;
d_front = z - dom.z_min;
d_back  = dom.z_max - z;
d_dom   = min([d_left, d_right, d_front, d_back]);

value      = [v_bottom; v_top; v_target; d_dom];
isterminal = [1;        1;     0;        1];
direction  = [-1;      -1;    -1;       -1];

end


%% Local JHTDB REST helper
function result = JHUTDB_getData(authToken, dataset, var_original, timepoint_original, ...
    temporal_method_original, spatial_method_original, ...
    spatial_operator_original, points, varargin)

    import matlab.net.http.*

    hasOption = ~isempty(varargin);
    if hasOption
        option = varargin{1};
    else
        option = [];
    end

    numPoints = size(points, 1);
    pointCells = arrayfun(@(i) sprintf('%.8f\t%.8f\t%.8f', ...
        points(i,1), points(i,2), points(i,3)), ...
        1:numPoints, 'UniformOutput', false);
    points_str = strjoin(pointCells, newline);

    options_http = HTTPOptions('ConnectTimeout', 1000);
    request = RequestMessage('POST', [], points_str);

    functionname = 'GetVariable';

    if ~hasOption
        url = ['https://web.idies.jhu.edu/turbulence-svc/values?authToken=', authToken, ...
               '&dataset=', dataset, ...
               '&function=', functionname, ...
               '&var=', var_original, ...
               '&t=', num2str(timepoint_original), ...
               '&sint=', spatial_method_original, ...
               '&sop=', spatial_operator_original, ...
               '&tint=', temporal_method_original];
    else
        if numel(option) < 2
            error('Optional input must be [time_end, delta_t].');
        end

        url = ['https://web.idies.jhu.edu/turbulence-svc/values?authToken=', authToken, ...
               '&dataset=', dataset, ...
               '&function=', functionname, ...
               '&var=', var_original, ...
               '&t=', num2str(timepoint_original), ...
               '&sint=', spatial_method_original, ...
               '&sop=', spatial_operator_original, ...
               '&tint=', temporal_method_original, ...
               '&timepoint_end=', num2str(option(1)), ...
               '&delta_t=', num2str(option(2))];
    end

    response = request.send(url, options_http);
    result = response.Body.Data;

    if response.StatusCode ~= matlab.net.http.StatusCode.OK
        if isstruct(result) && isfield(result, 'description')
            error(['HTTP Error ', char(response.StatusCode), ':', newline, ...
                   strjoin(result.description, newline)]);
        else
            error(['HTTP Error ', char(response.StatusCode), '.']);
        end
    end
end
