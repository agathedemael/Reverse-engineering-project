%% Blade Guard Shape Optimizer
%
% Optimization objective:
%   realMixingIndex = exchangeIngress * meanVorticity * upwardBounceFlow
%
% Hard geometry constraints:
%   - max radius R <= 1.00, same as the reference guards.
%   - total guard height <= 1.00, same as the hemisphere envelope.
%   - no candidate extends below the flat food-contact plane or above the dome.
%   - manufacturability filters: minimum rib width, minimum port height,
%     bounded recess depth, rounded ports, minimum wall thickness.

clear; clc; close all;
rng(23);

%% User settings
quickRun  = true;   
makeSTL   = true;
showPlots = true;

if quickRun
    searchCfg.Nx = 30; searchCfg.Ny = 30; searchCfg.Nz = 34;
    searchCfg.nSteps = 90; searchCfg.nJacobi = 35; searchCfg.Np = 2400;
    nCandidates = 30; nRefine = 5;
else
    searchCfg.Nx = 38; searchCfg.Ny = 38; searchCfg.Nz = 44;
    searchCfg.nSteps = 150; searchCfg.nJacobi = 55; searchCfg.Np = 3800;
    nCandidates = 90; nRefine = 8;
end

refineCfg.Nx = 50; refineCfg.Ny = 50; refineCfg.Nz = 58;
refineCfg.nSteps = 280; refineCfg.nJacobi = 85; refineCfg.Np = 6000;

common.dt = 0.0075;
common.nu = 0.0030;
common.domain_zmin = -0.12;       
common.domain_zmax = 1.08;
common.guard_zmin = 0.00;         
common.guard_zmax = 1.00;        
common.max_guard_radius = 1.00;
common.wall_clearance = 0.050;
common.shell_t = 0.050;
common.min_wall_thickness = 0.050;
common.min_rib_arc = 0.15;
common.min_opening_height = 0.13;
common.blade_z = 0.16;            % blade is near the flat base
common.blade_r = 0.24;
common.omega = 76;
common.startFraction = 0.40;

fprintf('\nBlade guard shape optimizer\n');
fprintf('Flat base at z=%.2f, dome top at z=%.2f, max radius R<=%.2f\n', ...
    common.guard_zmin, common.guard_zmax, common.max_guard_radius);
fprintf('Objective: exchange ingress x vorticity x upward bounce flow\n\n');

%% Generate candidates
candidates = generate_candidate_set(nCandidates, common);
fprintf('Generated %d manufacturable candidate shapes.\n\n', numel(candidates));

%% Stage 1 search
searchResults = repmat(empty_result(), numel(candidates), 1);
for k = 1:numel(candidates)
    fprintf('Stage 1 candidate %3d / %3d: %s\n', k, numel(candidates), candidates(k).name);
    searchResults(k) = simulate_candidate(candidates(k), searchCfg, common, false);
    print_result_line(searchResults(k));
end

[~, order] = sort([searchResults.realMixingIndex], 'descend');
topIdx = order(1:min(nRefine, numel(order)));

fprintf('\nStage 1 top candidates:\n');
for j = 1:numel(topIdx)
    r = searchResults(topIdx(j));
    fprintf('  %d. %s | realMix %.6g | exchange %.2f%% | vort %.4f | upward-bounce %.2f%%\n', ...
        j, r.name, r.realMixingIndex, 100*r.avgExchange, r.avgVort, 100*r.avgBounce);
end

%% Stage 2 refinement
refineResults = repmat(empty_result(), numel(topIdx), 1);
fprintf('\nRefining best %d candidates at higher resolution...\n\n', numel(topIdx));
for j = 1:numel(topIdx)
    k = topIdx(j);
    fprintf('Stage 2 refine %2d / %2d: %s\n', j, numel(topIdx), candidates(k).name);
    refineResults(j) = simulate_candidate(candidates(k), refineCfg, common, true);
    print_result_line(refineResults(j));
end

[~, bestLocal] = max([refineResults.realMixingIndex]);
bestResult = refineResults(bestLocal);
bestCandidate = bestResult.params;

%% Report
fprintf('\n================ OPTIMIZED CORRECT-ORIENTATION SHAPE ================\n');
fprintf('Best candidate:           %s\n', bestResult.name);
fprintf('Real mixing index:        %.8g\n', bestResult.realMixingIndex);
fprintf('Exchange ingress:         %.2f %%\n', 100*bestResult.avgExchange);
fprintf('Mean vorticity:           %.6f\n', bestResult.avgVort);
fprintf('Upward bounce flow:       %.2f %%\n', 100*bestResult.avgBounce);
fprintf('Top/dome reach:           %.2f %%\n', 100*bestResult.avgTopReach);
fprintf('Mean speed:               %.6f\n', bestResult.avgSpeed);
fprintf('\nOptimized geometry parameters:\n');
print_candidate(bestCandidate);
fprintf('=====================================================================\n\n');

Tsearch = result_table(searchResults);
Trefine = result_table(refineResults);
Tsearch = sortrows(Tsearch, 'realMixingIndex', 'descend');
Trefine = sortrows(Trefine, 'realMixingIndex', 'descend');

disp('Stage 1 ranking table:');
disp(Tsearch);
disp('Stage 2 refined ranking table:');
disp(Trefine);

%% Plotting results
if showPlots
    plot_best_shape(bestCandidate, common);
    plot_ranking(Trefine);
    plot_final_flow(bestResult, bestCandidate, common);
end

if makeSTL
    export_candidate_stl(bestCandidate, common, 'best_guard_optimized_correct_orientation.stl');
    fprintf('Exported STL: best_guard_optimized_correct_orientation.stl\n');
end

%% Local functions

function candidates = generate_candidate_set(nCandidates, common)
    candidates = struct([]);
    tries = 0;
    maxTries = 80*nCandidates;

    % Seed: broad lower rounded exchange ports plus shallow dome recesses.
    seed = candidate_from_values(7, 0.10, 0.30, 0.28, 0.12, 0.90, 0.0, 1.10, 0.11, 0.70, common);
    candidates = append_if_valid(candidates, seed, common);

    while numel(candidates) < nCandidates && tries < maxTries
        tries = tries + 1;
        nOpen = randi([5, 10]);
        spacing = 2*pi/nOpen;
        maxHalfWidth = min(0.40, 0.5*(spacing - common.min_rib_arc));
        minHalfWidth = 0.18;
        if maxHalfWidth <= minHalfWidth, continue; end

        recessDepth = 0.03 + 0.13*rand();
        openingHalfWidth = minHalfWidth + (maxHalfWidth-minHalfWidth)*rand();
        openingZCenter = 0.16 + 0.34*rand();       % lower-to-mid dome, not top
        openingHalfHeight = 0.070 + 0.080*rand();
        lipRadius = 0.86 + 0.10*rand();
        phase = 2*pi*rand();
        forcingGain = 0.85 + 0.45*rand();
        openingRadialHalfThickness = 0.090 + 0.055*rand();
        baseMouthBias = 0.45 + 0.40*rand();

        cand = candidate_from_values(nOpen, recessDepth, openingHalfWidth, ...
            openingZCenter, openingHalfHeight, lipRadius, phase, forcingGain, ...
            openingRadialHalfThickness, baseMouthBias, common);
        candidates = append_if_valid(candidates, cand, common);
    end

    if numel(candidates) < nCandidates
        warning('Only generated %d valid candidates.', numel(candidates));
    end

    for i = 1:numel(candidates)
        candidates(i).name = sprintf('Correct_%03d_N%d_D%.3f_W%.3f_Z%.2f', ...
            i, candidates(i).nOpenings, candidates(i).recessDepth, ...
            candidates(i).openingHalfWidth, candidates(i).openingZCenter);
    end
end

function cand = candidate_from_values(nOpen, recessDepth, openingHalfWidth, openingZCenter, openingHalfHeight, lipRadius, phase, forcingGain, openingRadialHalfThickness, baseMouthBias, common)
    cand = struct();
    cand.name = '';
    cand.nOpenings = nOpen;
    cand.recessDepth = recessDepth;
    cand.recessPower = 0.95;
    cand.openingHalfWidth = openingHalfWidth;
    cand.openingZCenter = openingZCenter;
    cand.openingHalfHeight = openingHalfHeight;
    cand.lipRadius = min(lipRadius, common.max_guard_radius);
    cand.phase = phase;
    cand.forcingGain = forcingGain;
    cand.openingRadialHalfThickness = openingRadialHalfThickness;
    cand.baseMouthBias = baseMouthBias;
    cand.wallThickness = common.shell_t;
    cand.maxRadius = common.max_guard_radius;
    cand.zMin = common.guard_zmin;
    cand.zMax = common.guard_zmax;
    cand.manufacturabilityScore = 0;
end

function candidates = append_if_valid(candidates, cand, common)
    [ok, score] = is_manufacturable(cand, common);
    if ok
        cand.manufacturabilityScore = score;
        candidates = [candidates; cand]; 
    end
end

function [ok, score] = is_manufacturable(c, common)
    spacing = 2*pi/c.nOpenings;
    ribArc = spacing - 2*c.openingHalfWidth;
    enoughRib = ribArc >= common.min_rib_arc;
    enoughOpeningHeight = 2*c.openingHalfHeight >= common.min_opening_height;
    boundedRadius = c.maxRadius <= common.max_guard_radius + 1e-12 && c.lipRadius <= common.max_guard_radius + 1e-12;
    boundedLength = c.zMin >= common.guard_zmin - 1e-12 && c.zMax <= common.guard_zmax + 1e-12;
    boundedRecess = c.recessDepth >= 0.00 && c.recessDepth <= 0.18;
    enoughWall = c.wallThickness >= common.min_wall_thickness - 1e-12;
    portInsideDome = (c.openingZCenter - c.openingHalfHeight) >= common.guard_zmin + 0.035 && ...
                     (c.openingZCenter + c.openingHalfHeight) <= 0.68;
    ok = enoughRib && enoughOpeningHeight && boundedRadius && boundedLength && ...
         boundedRecess && enoughWall && portInsideDome;

    ribScore = min(1, ribArc/0.28);
    wallScore = min(1, c.wallThickness/common.min_wall_thickness);
    recessScore = 1 - max(0, c.recessDepth - 0.14)/0.04;
    openingScore = min(1, (2*c.openingHalfHeight)/0.24);
    score = max(0, 0.35*ribScore + 0.25*wallScore + 0.20*recessScore + 0.20*openingScore);
end

function result = simulate_candidate(cand, cfg, common, keepFields)
    [X, Y, Z, x, y, z, dx, dy, dz, h, R] = make_grid(cfg, common);
    [fluid, solid, ~, exchangeReservoir] = build_domain_mask_param(cand, X, Y, Z, common);

    U = zeros(size(X)); V = zeros(size(X)); W = zeros(size(X)); P = zeros(size(X));
    [px, py, pz, ptag] = seed_particles_param(cfg.Np, fluid, exchangeReservoir, X, Y, Z, common);

    analysisRegion = fluid & Z > 0.02 & Z < 0.95 & R < 1.02;
    bounceRegion = fluid & R > 0.45 & R < 0.95 & Z > 0.18 & Z < 0.86;

    meanSpeed_t = zeros(cfg.nSteps,1);
    meanVort_t = zeros(cfg.nSteps,1);
    bounce_t = zeros(cfg.nSteps,1);
    topReach_t = zeros(cfg.nSteps,1);
    exchangeIn_t = zeros(cfg.nSteps,1);

    for n = 1:cfg.nSteps
        [Fx, Fy, Fz] = blade_forcing_param(X, Y, Z, common.blade_z, common.blade_r, common.omega, cand);
        U = U + common.dt * Fx;
        V = V + common.dt * Fy;
        W = W + common.dt * Fz;

        Uold = U; Vold = V; Wold = W;
        U = advect_field(Uold, Uold, Vold, Wold, X, Y, Z, common.dt, fluid);
        V = advect_field(Vold, Uold, Vold, Wold, X, Y, Z, common.dt, fluid);
        W = advect_field(Wold, Uold, Vold, Wold, X, Y, Z, common.dt, fluid);

        U = U + common.nu * common.dt * lap3d(U, dx, dy, dz);
        V = V + common.nu * common.dt * lap3d(V, dx, dy, dz);
        W = W + common.nu * common.dt * lap3d(W, dx, dy, dz);

        [U, V, W] = apply_noslip(U, V, W, solid);
        [U, V, W, P] = pressure_project(U, V, W, P, fluid, h, common.dt, cfg.nJacobi);
        [U, V, W] = apply_noslip(U, V, W, solid);

        domeCap = fluid & Z > 0.94;
        U(domeCap) = 0.985 * U(domeCap);
        V(domeCap) = 0.985 * V(domeCap);
        W(domeCap) = 0.985 * W(domeCap);

        [px, py, pz] = advect_particles(px, py, pz, U, V, W, X, Y, Z, common.dt, fluid);

        speed = sqrt(U.^2 + V.^2 + W.^2);
        [~, ~, ~, ommag] = curl3d(U, V, W, dx, dy, dz);

        meanSpeed_t(n) = safe_mean(speed(analysisRegion));
        meanVort_t(n) = safe_mean(ommag(analysisRegion));

        % Upward-bounce flow: positive vertical motion in the side/dome region.
        bounce_t(n) = safe_mean(W(bounceRegion) > 0);

        particleR = sqrt(px.^2 + py.^2);
        topReach_t(n) = mean(particleR > 0.38 & particleR < 0.90 & pz > 0.55 & pz < 0.98);

        exchangeOrigin = ptag == 2;
        if any(exchangeOrigin)
            exchangeIn_t(n) = mean(exchangeOrigin & particleR < 0.92 & pz > 0.00 & pz < 0.98);
        else
            exchangeIn_t(n) = 0;
        end
    end

    startIndex = max(1, round(common.startFraction * cfg.nSteps));
    result = empty_result();
    result.name = cand.name;
    result.params = cand;
    result.avgSpeed = mean(meanSpeed_t(startIndex:end));
    result.avgVort = mean(meanVort_t(startIndex:end));
    result.avgBounce = mean(bounce_t(startIndex:end));
    result.avgTopReach = mean(topReach_t(startIndex:end));
    result.avgExchange = mean(exchangeIn_t(startIndex:end));
    result.realMixingIndex = result.avgExchange * result.avgVort * result.avgBounce;
    result.manufacturabilityScore = cand.manufacturabilityScore;

    if keepFields
        result.finalU = U; result.finalV = V; result.finalW = W;
        result.finalFluid = fluid; result.finalSolid = solid;
        result.gridX = X; result.gridY = Y; result.gridZ = Z;
        result.x = x; result.y = y; result.z = z;
        result.px = px; result.py = py; result.pz = pz; result.ptag = ptag;
    end
end

function r = empty_result()
    r = struct('name','', 'params',struct(), 'avgSpeed',0, 'avgVort',0, ...
        'avgBounce',0, 'avgTopReach',0, 'avgExchange',0, 'realMixingIndex',0, ...
        'manufacturabilityScore',0, 'finalU',[], 'finalV',[], 'finalW',[], ...
        'finalFluid',[], 'finalSolid',[], 'gridX',[], 'gridY',[], 'gridZ',[], ...
        'x',[], 'y',[], 'z',[], 'px',[], 'py',[], 'pz',[], 'ptag',[]);
end

function T = result_table(results)
    n = numel(results);
    name = strings(n,1); realMixingIndex = zeros(n,1); avgExchange = zeros(n,1);
    avgVort = zeros(n,1); avgBounce = zeros(n,1); avgTopReach = zeros(n,1);
    avgSpeed = zeros(n,1); manufacturabilityScore = zeros(n,1);
    for i = 1:n
        name(i) = string(results(i).name);
        realMixingIndex(i) = results(i).realMixingIndex;
        avgExchange(i) = results(i).avgExchange;
        avgVort(i) = results(i).avgVort;
        avgBounce(i) = results(i).avgBounce;
        avgTopReach(i) = results(i).avgTopReach;
        avgSpeed(i) = results(i).avgSpeed;
        manufacturabilityScore(i) = results(i).manufacturabilityScore;
    end
    T = table(name, realMixingIndex, avgExchange, avgVort, avgBounce, avgTopReach, avgSpeed, manufacturabilityScore);
end

function [X, Y, Z, x, y, z, dx, dy, dz, h, R] = make_grid(cfg, common)
    x = linspace(-1.28, 1.28, cfg.Nx);
    y = linspace(-1.28, 1.28, cfg.Ny);
    z = linspace(common.domain_zmin, common.domain_zmax, cfg.Nz);
    [X, Y, Z] = ndgrid(x, y, z);
    dx = x(2)-x(1); dy = y(2)-y(1); dz = z(2)-z(1);
    h = mean([dx dy dz]);
    R = sqrt(X.^2 + Y.^2);
end

function [fluid, solid, openingFluid, exchangeReservoir] = build_domain_mask_param(c, X, Y, Z, common)
    R = sqrt(X.^2 + Y.^2);
    TH = atan2(Y, X);
    rWall = guard_radius_param(c, TH, Z, common);
    rFluidLimit = max(rWall - common.wall_clearance, 0);

    insideDome = Z >= common.guard_zmin & Z <= common.guard_zmax & R < rFluidLimit;
    foodReservoir = Z >= common.domain_zmin & Z < common.guard_zmin & R < c.lipRadius - common.wall_clearance;

    openingFluid = false(size(X));
    centers = c.phase + (0:c.nOpenings-1)*(2*pi/c.nOpenings);
    for hc = centers
        dth = atan2(sin(TH - hc), cos(TH - hc));
        zNorm = (Z - c.openingZCenter) / c.openingHalfHeight;
        thNorm = dth / c.openingHalfWidth;
        roundedPort = (zNorm.^2 + thNorm.^2) < 1.0;
        nearWall = abs(R - rWall) < c.openingRadialHalfThickness;
        openingFluid = openingFluid | (roundedPort & nearWall);

        % Smooth continuation into the flat base, so ports are printable and
        % can draw food from the contact plane without creating sharp slots.
        baseMouth = abs(dth) < c.baseMouthBias*c.openingHalfWidth & ...
                    Z >= common.guard_zmin & Z < min(0.16, c.openingZCenter + c.openingHalfHeight) & ...
                    R > 0.68 & R < 1.00;
        openingFluid = openingFluid | baseMouth;
    end

    exchangeReservoir = (Z > common.domain_zmin + 0.02 & Z < 0.28 & R > 0.78 & R < 1.15) | ...
                        (Z >= common.guard_zmin & Z < 0.55 & R > max(rWall + common.wall_clearance, 0.74) & R < 1.15);

    fluid = insideDome | foodReservoir | openingFluid | exchangeReservoir;
    fluid = fluid & Z >= common.domain_zmin & Z <= common.domain_zmax;

    domeWall = Z >= common.guard_zmin & Z <= common.guard_zmax & abs(R - rWall) <= common.shell_t;
    domeWall(openingFluid) = false;
    baseContactRing = Z >= -0.01 & Z <= 0.025 & R > 0.82 & R < 1.00;
    baseContactRing(openingFluid) = false;
    outerReservoirWall = Z > common.domain_zmin+0.02 & Z < 0.55 & abs(R - 1.15) <= common.shell_t;
    shaft = R < 0.055 & Z > 0.10 & Z < common.guard_zmax;

    solid = ~fluid | domeWall | baseContactRing | outerReservoirWall | shaft;
    solid(1,:,:) = true; solid(end,:,:) = true;
    solid(:,1,:) = true; solid(:,end,:) = true;
    solid(:,:,1) = true; solid(:,:,end) = true;
    fluid(solid) = false;
    openingFluid = openingFluid & fluid;
    exchangeReservoir = exchangeReservoir & fluid;
end

function rWall = guard_radius_param(c, TH, Z, common)
    % Hemisphere with flat/open base at z=0 and rounded top at z=1:
    % r(z) = sqrt(1 - z^2). This has max radius at the flat food side.
    zeta = min(max((Z-common.guard_zmin)/(common.guard_zmax-common.guard_zmin), 0), 1);
    baseWallR = sqrt(max(1 - zeta.^2, 0));

    % Recesses are strongest near lower/mid dome and fade near top
    fadeTop = (1 - zeta).^c.recessPower;
    lobe = 0.5 + 0.5*cos(c.nOpenings*(TH - c.phase));
    recess = c.recessDepth .* lobe .* fadeTop;
    rWall = baseWallR .* (1 - recess);
    rWall = min(rWall, common.max_guard_radius);
    rWall = max(rWall, 0);
end

function [Fx, Fy, Fz] = blade_forcing_param(X, Y, Z, blade_z, blade_r, omega, c)
    R = sqrt(X.^2 + Y.^2) + 1e-8;
    TH = atan2(Y, X);

    axialCore = exp(-((Z - blade_z) / 0.080).^2);
    Ftheta = omega * (R / blade_r) .* exp(-(R / (blade_r * 1.15)).^2) .* axialCore;

    % The artificial blade pulls material from the flat base
    % and throws it upward into the dome.
    Fz_core_up = 18.0 .* exp(-((Z - (blade_z + 0.02)) / 0.15).^2) .* exp(-(R / (blade_r * 1.35)).^2);
    Fr_lower_out = 8.0 .* exp(-((Z - 0.18) / 0.16).^2) .* exp(-((R - blade_r) / 0.22).^2);

    % Near the rounded dome/edge, the flow bounces/redirects: upward motion
    % near the side wall plus inward curling near the upper curvature.
    sideUp = exp(-((R - 0.68) / 0.20).^2) .* exp(-((Z - 0.48) / 0.34).^2);
    Fz_side_up = 13.0 * sideUp;
    upperCurve = exp(-((R - 0.50) / 0.22).^2) .* exp(-((Z - 0.78) / 0.16).^2);
    Fr_upper_in = -9.0 * upperCurve;
    Fz_upper_slow = -3.0 * upperCurve; 

    centerDownReturn = exp(-(R / 0.34).^2) .* exp(-((Z - 0.70) / 0.22).^2);
    Fz_center_down = -5.0 * centerDownReturn;

    Fr = Fr_lower_out + Fr_upper_in;
    Fz = Fz_core_up + Fz_side_up + Fz_upper_slow + Fz_center_down;

    centers = c.phase + (0:c.nOpenings-1)*(2*pi/c.nOpenings);
    exchangeZone = zeros(size(X));
    for hc = centers
        dth = atan2(sin(TH - hc), cos(TH - hc));
        exchangeZone = exchangeZone + exp(-(dth / max(c.openingHalfWidth,0.08)).^2) .* ...
            exp(-((Z - c.openingZCenter) / max(c.openingHalfHeight,0.06)).^2) .* ...
            exp(-((R - 0.93) / 0.18).^2);
    end
    Fr_exchange_in = -7.0 * c.forcingGain * exchangeZone;
    Fz_exchange_up = 6.0 * c.forcingGain * exchangeZone;
    swirlLift = 2.5 * c.forcingGain * exchangeZone .* sin(c.nOpenings*(TH - c.phase));

    Fr = Fr + Fr_exchange_in;
    Fz = Fz + Fz_exchange_up + swirlLift;

    Fx = Fr .* cos(TH) - Ftheta .* sin(TH);
    Fy = Fr .* sin(TH) + Ftheta .* cos(TH);
end

function phi_new = advect_field(phi, U, V, W, X, Y, Z, dt, fluid)
    Xb = X - dt * U; Yb = Y - dt * V; Zb = Z - dt * W;
    phi_new = interpn(X, Y, Z, phi, Xb, Yb, Zb, 'linear', 0);
    phi_new(~fluid) = 0;
end

function L = lap3d(A, dx, dy, dz)
    L = (circshift(A,[-1 0 0])-2*A+circshift(A,[1 0 0]))/dx^2 + ...
        (circshift(A,[0 -1 0])-2*A+circshift(A,[0 1 0]))/dy^2 + ...
        (circshift(A,[0 0 -1])-2*A+circshift(A,[0 0 1]))/dz^2;
end

function [U, V, W] = apply_noslip(U, V, W, solid)
    U(solid)=0; V(solid)=0; W(solid)=0;
    U(1,:,:)=0; U(end,:,:)=0; U(:,1,:)=0; U(:,end,:)=0; U(:,:,1)=0; U(:,:,end)=0;
    V(1,:,:)=0; V(end,:,:)=0; V(:,1,:)=0; V(:,end,:)=0; V(:,:,1)=0; V(:,:,end)=0;
    W(1,:,:)=0; W(end,:,:)=0; W(:,1,:)=0; W(:,end,:)=0; W(:,:,1)=0; W(:,:,end)=0;
end

function [U, V, W, P] = pressure_project(U, V, W, P, fluid, h, dt, nIter)
    div = divergence3d(U, V, W, h, h, h);
    div(~fluid) = 0;
    for it = 1:nIter
        Pe = circshift(P,[-1 0 0]); Pw = circshift(P,[1 0 0]);
        Pn = circshift(P,[0 -1 0]); Ps = circshift(P,[0 1 0]);
        Pt = circshift(P,[0 0 -1]); Pb = circshift(P,[0 0 1]);
        Fe = circshift(fluid,[-1 0 0]); Fw = circshift(fluid,[1 0 0]);
        Fn = circshift(fluid,[0 -1 0]); Fs = circshift(fluid,[0 1 0]);
        Ft = circshift(fluid,[0 0 -1]); Fb = circshift(fluid,[0 0 1]);
        sumP = Pe.*Fe + Pw.*Fw + Pn.*Fn + Ps.*Fs + Pt.*Ft + Pb.*Fb;
        cnt = double(Fe)+double(Fw)+double(Fn)+double(Fs)+double(Ft)+double(Fb);
        cnt(cnt==0)=1;
        P = (sumP - h^2*(div/dt))./cnt;
        P(~fluid)=0;
    end
    [dPdx, dPdy, dPdz] = grad3d(P, h, h, h);
    U = U - dt*dPdx; V = V - dt*dPdy; W = W - dt*dPdz;
    U(~fluid)=0; V(~fluid)=0; W(~fluid)=0;
end

function div = divergence3d(U, V, W, dx, dy, dz)
    dUdx = (circshift(U,[-1 0 0])-circshift(U,[1 0 0]))/(2*dx);
    dVdy = (circshift(V,[0 -1 0])-circshift(V,[0 1 0]))/(2*dy);
    dWdz = (circshift(W,[0 0 -1])-circshift(W,[0 0 1]))/(2*dz);
    div = dUdx + dVdy + dWdz;
end

function [dAdx, dAdy, dAdz] = grad3d(A, dx, dy, dz)
    dAdx = (circshift(A,[-1 0 0])-circshift(A,[1 0 0]))/(2*dx);
    dAdy = (circshift(A,[0 -1 0])-circshift(A,[0 1 0]))/(2*dy);
    dAdz = (circshift(A,[0 0 -1])-circshift(A,[0 0 1]))/(2*dz);
end

function [omx, omy, omz, ommag] = curl3d(U, V, W, dx, dy, dz)
    [~, dUdy, dUdz] = grad3d(U, dx, dy, dz);
    [dVdx, ~, dVdz] = grad3d(V, dx, dy, dz);
    [dWdx, dWdy, ~] = grad3d(W, dx, dy, dz);
    omx = dWdy - dVdz;
    omy = dUdz - dWdx;
    omz = dVdx - dUdy;
    ommag = sqrt(omx.^2 + omy.^2 + omz.^2);
end

function [px, py, pz, ptag] = seed_particles_param(Np, fluid, exchangeReservoir, X, Y, Z, common)
    R = sqrt(X.^2 + Y.^2);
    foodBase = fluid & Z > -0.10 & Z < 0.10 & R < 0.88;
    bladeRegion = fluid & Z > 0.08 & Z < 0.32 & R < 0.72;
    domeRegion = fluid & Z > 0.45 & Z < 0.95 & R > 0.20 & R < 0.88;
    exchangeFood = exchangeReservoir;

    nExchange = round(0.30*Np);
    nBase = round(0.34*Np);
    nBlade = round(0.24*Np);
    nDome = Np - nExchange - nBase - nBlade;

    idxBase = find(foodBase); if isempty(idxBase), idxBase = find(fluid); end
    idxBlade = find(bladeRegion); if isempty(idxBlade), idxBlade = find(fluid); end
    idxDome = find(domeRegion); if isempty(idxDome), idxDome = find(fluid); end
    idxExchange = find(exchangeFood); if isempty(idxExchange), idxExchange = find(fluid); end

    pickBase = idxBase(randi(numel(idxBase), [nBase 1]));
    pickBlade = idxBlade(randi(numel(idxBlade), [nBlade 1]));
    pickDome = idxDome(randi(numel(idxDome), [nDome 1]));
    pickExchange = idxExchange(randi(numel(idxExchange), [nExchange 1]));

    pick = [pickBase; pickBlade; pickDome; pickExchange];
    ptag = [ones(nBase+nBlade+nDome,1); 2*ones(nExchange,1)];
    px = X(pick); py = Y(pick); pz = Z(pick);
end

function [px, py, pz] = advect_particles(px, py, pz, U, V, W, X, Y, Z, dt, fluid)
    Ui = interpn(X,Y,Z,U,px,py,pz,'linear',0);
    Vi = interpn(X,Y,Z,V,px,py,pz,'linear',0);
    Wi = interpn(X,Y,Z,W,px,py,pz,'linear',0);
    turb = 0.017;
    pxt = px + dt*Ui + turb*randn(size(px))*sqrt(dt);
    pyt = py + dt*Vi + turb*randn(size(py))*sqrt(dt);
    pzt = pz + dt*Wi + 0.70*turb*randn(size(pz))*sqrt(dt);
    valid = interpn(X,Y,Z,double(fluid),pxt,pyt,pzt,'nearest',0) > 0.5;
    px(valid)=pxt(valid); py(valid)=pyt(valid); pz(valid)=pzt(valid);
    hit = ~valid;
    if any(hit)
        px(hit) = px(hit) - 0.40*dt*Ui(hit) + 0.006*randn(sum(hit),1);
        py(hit) = py(hit) - 0.40*dt*Vi(hit) + 0.006*randn(sum(hit),1);
        pz(hit) = pz(hit) - 0.40*dt*Wi(hit) + 0.004*randn(sum(hit),1);
    end
    stillValid = interpn(X,Y,Z,double(fluid),px,py,pz,'nearest',0) > 0.5;
    if any(~stillValid)
        idxFluid = find(fluid);
        bad = find(~stillValid);
        replacement = idxFluid(randi(numel(idxFluid), [numel(bad) 1]));
        px(bad)=X(replacement); py(bad)=Y(replacement); pz(bad)=Z(replacement);
    end
end

function m = safe_mean(A)
    A = A(:);
    if isempty(A), m = 0; else, m = mean(A); end
end

function print_result_line(r)
    fprintf('  realMix %.6g | exchange %.2f%% | vort %.4f | upward-bounce %.2f%% | manuf %.2f\n', ...
        r.realMixingIndex, 100*r.avgExchange, r.avgVort, 100*r.avgBounce, r.manufacturabilityScore);
end

function print_candidate(c)
    fprintf('  nOpenings:                  %d\n', c.nOpenings);
    fprintf('  inward recess depth:         %.4f of local radius\n', c.recessDepth);
    fprintf('  opening half-width:          %.4f rad  full %.4f rad\n', c.openingHalfWidth, 2*c.openingHalfWidth);
    fprintf('  opening z-center:            %.4f  [0=flat base, 1=dome top]\n', c.openingZCenter);
    fprintf('  opening full height:         %.4f\n', 2*c.openingHalfHeight);
    fprintf('  lip/base radius:             %.4f\n', c.lipRadius);
    fprintf('  wall thickness:              %.4f\n', c.wallThickness);
    fprintf('  phase angle:                 %.4f rad\n', c.phase);
    fprintf('  exchange forcing gain:       %.4f\n', c.forcingGain);
    fprintf('  radial opening half-thick.:  %.4f\n', c.openingRadialHalfThickness);
    fprintf('  base-mouth bias:             %.4f\n', c.baseMouthBias);
    fprintf('  manufacturability score:     %.3f / 1.000\n', c.manufacturabilityScore);
end

function plot_best_shape(c, common)
    [Xs, Ys, Zs] = candidate_surface(c, common, 150);
    figure('Color','w','Name','Optimized correct-orientation blade guard','Position',[140 80 1050 820]);
    surf(Xs,Ys,Zs,'FaceAlpha',0.30,'EdgeColor','none');
    axis equal; grid on; xlabel('x'); ylabel('y'); zlabel('z');
    title('Optimized guard: flat base down, rounded dome up');
    view(38,24); camlight; lighting gouraud; hold on;
    th = linspace(0,2*pi,240);
    plot3(cos(th), sin(th), zeros(size(th)), 'k--', 'LineWidth', 1.2);
    text(0,0,-0.055,'flat food-contact base, R=1 envelope','HorizontalAlignment','center');
end

function [Xs, Ys, Zs] = candidate_surface(c, common, res)
    [TH, ZZ] = meshgrid(linspace(0,2*pi,res), linspace(common.guard_zmin, common.guard_zmax,res));
    zeta = ZZ;
    baseR = sqrt(max(1-zeta.^2,0));
    fadeTop = (1-zeta).^c.recessPower;
    lobe = 0.5 + 0.5*cos(c.nOpenings*(TH-c.phase));
    Rs = baseR .* (1 - c.recessDepth .* lobe .* fadeTop);
    Xs = Rs.*cos(TH); Ys = Rs.*sin(TH); Zs = ZZ;

    mask = true(size(Xs));
    centers = c.phase + (0:c.nOpenings-1)*(2*pi/c.nOpenings);
    for hc = centers
        dth = atan2(sin(TH-hc), cos(TH-hc));
        zNorm = (Zs - c.openingZCenter)/c.openingHalfHeight;
        thNorm = dth/c.openingHalfWidth;
        sideHole = (zNorm.^2 + thNorm.^2) < 1.0;
        baseHole = abs(dth) < c.baseMouthBias*c.openingHalfWidth & Zs < 0.16 & Rs > 0.68;
        mask = mask & ~(sideHole | baseHole);
    end
    Xs(~mask)=nan; Ys(~mask)=nan; Zs(~mask)=nan;
end

function plot_ranking(T)
    figure('Color','w','Name','Refined optimization ranking','Position',[220 180 1000 520]);
    bar(T.realMixingIndex);
    ylabel('Real mixing efficiency index');
    title('Refined candidates ranked by corrected-orientation real mixing index');
    xticks(1:height(T)); xticklabels(T.name); xtickangle(25); grid on;
end

function plot_final_flow(result, cand, common)
    if isempty(result.finalU), return; end
    X = result.gridX; Y = result.gridY; Z = result.gridZ;
    x = result.x; y = result.y; z = result.z;
    U = result.finalU; V = result.finalV; W = result.finalW;
    speed = sqrt(U.^2 + V.^2 + W.^2);
    [~,iy0] = min(abs(y));
    [~,izBlade] = min(abs(z - common.blade_z));

    figure('Color','w','Name','Correct-orientation optimized CFD field','Position',[80 80 1400 760]);
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

    nexttile;
    [Xs,Ys,Zs] = candidate_surface(cand, common, 120);
    surf(Xs,Ys,Zs,'FaceAlpha',0.16,'EdgeColor','none'); hold on;
    sub = randperm(numel(result.px), min(2500, numel(result.px)));
    scatter3(result.px(sub), result.py(sub), result.pz(sub), 10, result.ptag(sub), 'filled');
    axis equal; grid on; view(35,22); xlabel('x'); ylabel('y'); zlabel('z');
    title('Optimized guard and tracer particles'); camlight; lighting gouraud;

    nexttile;
    imagesc(x,z,squeeze(speed(:,iy0,:))'); set(gca,'YDir','normal'); axis equal tight;
    colorbar; title('Speed magnitude, XZ center slice'); xlabel('x'); ylabel('z'); hold on;
    [bx,bz] = candidate_outline_xz(cand, common, 350);
    plot(bx,bz,'k','LineWidth',1.6); plot(-bx,bz,'k','LineWidth',1.6);
    plot([-1 1],[0 0],'k--','LineWidth',1.1);

    nexttile;
    imagesc(x,y,squeeze(speed(:,:,izBlade))'); set(gca,'YDir','normal'); axis equal tight;
    colorbar; title('Speed magnitude near blade / flat base'); xlabel('x'); ylabel('y'); hold on;
    th = linspace(0,2*pi,220);
    plot(common.blade_r*cos(th),common.blade_r*sin(th),'k--','LineWidth',1.2);
end

function [xline, zline] = candidate_outline_xz(c, common, n)
    zline = linspace(common.guard_zmin, common.guard_zmax, n);
    zeta = zline;
    baseR = sqrt(max(1-zeta.^2,0));
    fadeTop = (1-zeta).^c.recessPower;
    lobe = 0.5 + 0.5*cos(c.nOpenings*(0-c.phase));
    xline = baseR .* (1 - c.recessDepth .* lobe .* fadeTop);
end

function export_candidate_stl(c, common, filename)
    [X,Y,Z] = candidate_surface(c, common, 130);
    [uu,vv] = meshgrid(linspace(0,1,size(X,2)), linspace(0,1,size(X,1)));
    valid = ~isnan(X);
    UV = [uu(valid), vv(valid)];
    P = [X(valid), Y(valid), Z(valid)];
    [UVu, ia, ~] = unique(UV,'rows');
    P = P(ia,:);
    DT = delaunayTriangulation(UVu);
    TR = triangulation(DT.ConnectivityList, P);
    write_stl_ascii(filename, TR);
end

function write_stl_ascii(filename, TR)
    F = TR.ConnectivityList; V = TR.Points;
    fid = fopen(filename,'w');
    fprintf(fid,'solid optimized_blade_guard_correct_orientation\n');
    for i = 1:size(F,1)
        v1 = V(F(i,1),:); v2 = V(F(i,2),:); v3 = V(F(i,3),:);
        n = cross(v2-v1, v3-v1);
        if norm(n)>0, n = n/norm(n); else, n = [0 0 0]; end
        fprintf(fid,' facet normal %.6f %.6f %.6f\n', n);
        fprintf(fid,'  outer loop\n');
        fprintf(fid,'   vertex %.6f %.6f %.6f\n', v1);
        fprintf(fid,'   vertex %.6f %.6f %.6f\n', v2);
        fprintf(fid,'   vertex %.6f %.6f %.6f\n', v3);
        fprintf(fid,'  endloop\n');
        fprintf(fid,' endfacet\n');
    end
    fprintf(fid,'endsolid optimized_blade_guard_correct_orientation\n');
    fclose(fid);
end
