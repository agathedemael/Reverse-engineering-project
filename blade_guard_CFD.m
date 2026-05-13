clear; clc; close all; rng(3);

%% Grid
Nx = 48;
Ny = 48;
Nz = 56;
x = linspace(-1.28, 1.28, Nx);
y = linspace(-1.28, 1.28, Ny);
z = linspace(-1.36, -0.12, Nz);
[X, Y, Z] = ndgrid(x, y, z);
dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);
h  = mean([dx dy dz]);
R  = sqrt(X.^2 + Y.^2);
TH = atan2(Y, X);

%% Simulation settings
dt      = 0.0085;
nSteps  = 260;
nu      = 0.0032;
nJacobi = 80;

%% Blade settings
blade_z = -0.38;
blade_r = 0.23;
omega   = 72;

%% Hybrid blade guard domain settings
domain_zmin    = -1.35;
domain_zmax    = -0.20;
wall_clearance = 0.055;

model_names = { ...
    '1. Baseline Closed', ...
    '2. Perforated / Side Slits', ...
    '3. Daisy / Bottom Exchange'};

nModels = 3;
all_speed       = zeros(nSteps, nModels);
all_vort        = zeros(nSteps, nModels);
all_recirc      = zeros(nSteps, nModels);
all_topCurve    = zeros(nSteps, nModels);
all_exchangeIn  = zeros(nSteps, nModels);

%% Run each guard model

for model = 1:nModels

    fprintf('\nRunning model %d / %d: %s\n', model, nModels, model_names{model});
    [fluid, solid, bellSolid, openingFluid, exchangeReservoir] = ...
        build_domain_mask_hybrid(model, X, Y, Z, ...
        domain_zmin, domain_zmax, wall_clearance);

    U = zeros(size(X));
    V = zeros(size(X));
    W = zeros(size(X));
    P = zeros(size(X));
    Np = 4800;
    [px, py, pz, ptag] = seed_particles_hybrid(Np, fluid, exchangeReservoir, X, Y, Z, model);

    analysisRegion = fluid & ...
        Z > domain_zmin & ...
        Z < domain_zmax & ...
        R < 1.08;

    recircRegion = fluid & ...
        R > 0.42 & R < 0.92 & ...
        Z > -1.15 & Z < -0.35;

    topCurveRegion = fluid & ...
        R > 0.45 & R < 0.88 & ...
        Z > -0.45 & Z < -0.22;

    innerBellRegion = fluid & ...
        R < 0.88 & ...
        Z > -1.20 & Z < -0.22;

    meanSpeed_t   = zeros(nSteps, 1);
    meanVort_t    = zeros(nSteps, 1);
    recirc_t      = zeros(nSteps, 1);
    topCurve_t    = zeros(nSteps, 1);
    exchangeIn_t  = zeros(nSteps, 1);

    % Time integration
    for n = 1:nSteps

        %% Blade forcing 
        [Fx, Fy, Fz] = blade_forcing_hybrid(X, Y, Z, blade_z, blade_r, omega, model);

        U = U + dt * Fx;
        V = V + dt * Fy;
        W = W + dt * Fz;

        %% Semi-Lagrangian advection
        Uold = U;
        Vold = V;
        Wold = W;

        U = advect_field(Uold, Uold, Vold, Wold, X, Y, Z, dt, fluid);
        V = advect_field(Vold, Uold, Vold, Wold, X, Y, Z, dt, fluid);
        W = advect_field(Wold, Uold, Vold, Wold, X, Y, Z, dt, fluid);

        %% Viscosity / diffusion
        U = U + nu * dt * lap3d(U, dx, dy, dz);
        V = V + nu * dt * lap3d(V, dx, dy, dz);
        W = W + nu * dt * lap3d(W, dx, dy, dz);

        %% No-slip at solid walls
        [U, V, W] = apply_noslip(U, V, W, solid);

        %% Pressure projection
        [U, V, W, P] = pressure_project(U, V, W, P, fluid, h, dt, nJacobi);

        %% Re-apply no-slip
        [U, V, W] = apply_noslip(U, V, W, solid);

        %% Mild top damping
        topBand = fluid & Z > -0.245;
        U(topBand) = 0.990 * U(topBand);
        V(topBand) = 0.990 * V(topBand);
        W(topBand) = 0.990 * W(topBand);

        %% Particle advection with wall rebound
        [px, py, pz] = advect_particles(px, py, pz, U, V, W, X, Y, Z, dt, fluid);

        speed = sqrt(U.^2 + V.^2 + W.^2);
        [~, ~, ~, ommag] = curl3d(U, V, W, dx, dy, dz);

        meanSpeed_t(n) = safe_mean(speed(analysisRegion));
        meanVort_t(n)  = safe_mean(ommag(analysisRegion));
        recirc_t(n)    = safe_mean(W(recircRegion) > 0);

        particleTopCurve = ...
            sqrt(px.^2 + py.^2) > 0.45 & ...
            sqrt(px.^2 + py.^2) < 0.88 & ...
            pz > -0.45 & ...
            pz < -0.22;

        topCurve_t(n) = mean(particleTopCurve);
        exchangeOrigin = ptag == 2;

        if any(exchangeOrigin)
            exchangeIn_t(n) = mean( ...
                exchangeOrigin & ...
                sqrt(px.^2 + py.^2) < 0.88 & ...
                pz > -1.20 & ...
                pz < -0.22);
        else
            exchangeIn_t(n) = 0;
        end

        if mod(n, 20) == 0
            fprintf('  step %3d / %3d | top curve %.1f%% | exchange ingress %.1f%%\n', ...
                n, nSteps, 100*topCurve_t(n), 100*exchangeIn_t(n));
        end
    end

    all_speed(:, model)      = meanSpeed_t;
    all_vort(:, model)       = meanVort_t;
    all_recirc(:, model)     = recirc_t;
    all_topCurve(:, model)   = topCurve_t;
    all_exchangeIn(:, model) = exchangeIn_t;

    speed = sqrt(U.^2 + V.^2 + W.^2);
    [~, ~, ~, ommag] = curl3d(U, V, W, dx, dy, dz);
    iy0 = round(Ny / 2);
    [~, izBlade] = min(abs(z - blade_z));

    % Plots and figures
    figure('Color', 'w', 'Name', model_names{model}, 'Position', [60 50 1750 980]);
    tl = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    title(tl, ['Hybrid Blade-Guard CFD Flow Analysis — ' model_names{model}], ...
        'FontWeight', 'bold', 'FontSize', 14);

    %% Particles and blade guard
    nexttile;
    hold on;

    [Xs, Ys, Zs] = bell_surface(model, 120);

    surf(Xs, Ys, Zs, ...
        'FaceColor', [0.45 0.75 0.95], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.18);

    draw_inside_top_curve_region();

    if model == 2 || model == 3
        draw_exchange_reservoir_region(model);
    end

    sub = randperm(numel(px), min(2600, numel(px)));
    scatter3(px(sub), py(sub), pz(sub), 13, ptag(sub), 'filled');
    th = linspace(0, 2*pi, 220);
    plot3(blade_r*cos(th), blade_r*sin(th), blade_z*ones(size(th)), ...
        'k--', 'LineWidth', 1.5);

    xlabel('x');
    ylabel('y');
    zlabel('z');

    if model == 2
        title('3D Tracers: side food enters through slits');
    elseif model == 3
        title('3D Tracers: lower food enters through daisy rim');
    else
        title('3D Tracers: closed internal guard circulation');
    end

    axis equal;
    xlim([-1.18 1.18]);
    ylim([-1.18 1.18]);
    zlim([-1.35 -0.16]);
    view(35, 22);
    grid on;
    camlight;
    lighting gouraud;

    cb = colorbar;
    ylabel(cb, 'Particle source: 1 = inside, 2 = exchange reservoir');

    %% XZ slice
    nexttile;
    imagesc(x, z, squeeze(speed(:, iy0, :))');
    set(gca, 'YDir', 'normal');
    axis equal tight;
    colormap turbo;
    colorbar;
    xlabel('x');
    ylabel('z');
    title('Speed Magnitude: XZ Slice');
    hold on;

    [bx, bz] = bell_outline_xz(model, 350);
    plot(bx, bz, 'k', 'LineWidth', 1.7);
    plot(-bx, bz, 'k', 'LineWidth', 1.7);

    rectangle('Position', [-0.88 -0.45 1.76 0.23], ...
        'EdgeColor', 'w', ...
        'LineWidth', 1.4, ...
        'LineStyle', '--');

    if model == 2
        rectangle('Position', [-1.14 -1.04 2.28 0.30], ...
            'EdgeColor', 'c', ...
            'LineWidth', 1.4, ...
            'LineStyle', '--');

        text(-1.12, -0.72, 'side slit exchange zone', ...
            'Color', 'c', ...
            'FontWeight', 'bold', ...
            'FontSize', 9);

    elseif model == 3
        rectangle('Position', [-1.15 -1.30 2.30 0.36], ...
            'EdgeColor', 'c', ...
            'LineWidth', 1.4, ...
            'LineStyle', '--');

        text(-1.12, -0.94, 'bottom daisy-rim exchange zone', ...
            'Color', 'c', ...
            'FontWeight', 'bold', ...
            'FontSize', 9);
    end

    %% XY slice near the blade
    nexttile;
    imagesc(x, y, squeeze(speed(:, :, izBlade))');
    set(gca, 'YDir', 'normal');
    axis equal tight;
    colormap turbo;
    colorbar;
    xlabel('x');
    ylabel('y');
    title(sprintf('Vortex Speed: XY Slice at z = %.2f', z(izBlade)));
    hold on;
    plot(blade_r*cos(th), blade_r*sin(th), 'k--', 'LineWidth', 1.2);

    %% XZ vorticity slice
    nexttile;
    imagesc(x, z, squeeze(ommag(:, iy0, :))');
    set(gca, 'YDir', 'normal');
    axis equal tight;
    colormap turbo;
    colorbar;
    xlabel('x');
    ylabel('z');
    title('Vorticity Magnitude: XZ Slice');
    hold on;
    plot(bx, bz, 'w', 'LineWidth', 1.4);
    plot(-bx, bz, 'w', 'LineWidth', 1.4);

    %% Velocity vector field in XZ
    nexttile;
    xs = squeeze(X(:, iy0, :));
    zs = squeeze(Z(:, iy0, :));
    us = squeeze(U(:, iy0, :));
    ws = squeeze(W(:, iy0, :));
    ss = squeeze(speed(:, iy0, :));

    imagesc(x, z, ss');
    set(gca, 'YDir', 'normal');
    hold on;

    skip = 3;
    quiver(xs(1:skip:end, 1:skip:end), ...
           zs(1:skip:end, 1:skip:end), ...
           us(1:skip:end, 1:skip:end), ...
           ws(1:skip:end, 1:skip:end), ...
           1.8, 'k');

    axis equal tight;
    colormap turbo;
    colorbar;
    xlabel('x');
    ylabel('z');

    if model == 2
        title('Velocity Vectors: suction relief through side slits');
    elseif model == 3
        title('Velocity Vectors: bottom intake through daisy rim');
    else
        title('Velocity Vectors: closed internal circulation');
    end

    plot(bx, bz, 'w', 'LineWidth', 1.4);
    plot(-bx, bz, 'w', 'LineWidth', 1.4);

    text(-1.04, -0.86, 'upward wall return', ...
        'Color', 'w', 'FontWeight', 'bold', 'FontSize', 9);

    text(-0.70, -0.30, 'upper roll-over', ...
        'Color', 'w', 'FontWeight', 'bold', 'FontSize', 9);

    text(-0.36, -0.55, 'blade suction', ...
        'Color', 'w', 'FontWeight', 'bold', 'FontSize', 9);

    if model == 2
        text(0.62, -0.88, 'side slit inflow', ...
            'Color', 'c', 'FontWeight', 'bold', 'FontSize', 9);
    elseif model == 3
        text(0.48, -1.12, 'bottom rim intake', ...
            'Color', 'c', 'FontWeight', 'bold', 'FontSize', 9);
    end

    %% Time-history metrics
    nexttile;
    t = (1:nSteps) * dt;

    yyaxis left;
    plot(t, meanSpeed_t, 'LineWidth', 1.8);
    hold on;
    plot(t, meanVort_t, 'LineWidth', 1.8);
    ylabel('Mean speed / mean vorticity');

    yyaxis right;
    plot(t, 100 * topCurve_t, '--', 'LineWidth', 2.0);
    hold on;
    plot(t, 100 * exchangeIn_t, ':', 'LineWidth', 2.3);
    ylabel('Particle fraction (%)');

    xlabel('Time');
    title('Internal Flow Metrics');
    legend('Mean speed', ...
           'Mean vorticity', ...
           'Upper-curve particle fraction', ...
           'Exchange-origin particles entering blade guard', ...
           'Location', 'best');
    grid on;

end

% Comparison figure
t = (1:nSteps) * dt;
figure('Color', 'w', 'Position', [120 70 1280 950]);
tl2 = tiledlayout(5, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl2, 'Hybrid Blade Guard Flow Comparison', ...
    'FontWeight', 'bold', 'FontSize', 14);

nexttile;
plot(t, all_speed(:,1), 'LineWidth', 2);
hold on;
plot(t, all_speed(:,2), 'LineWidth', 2);
plot(t, all_speed(:,3), 'LineWidth', 2);
grid on;
ylabel('Mean speed');
title('Mean Speed in Solved Food Domain');
legend(model_names, 'Location', 'best');

nexttile;
plot(t, all_vort(:,1), 'LineWidth', 2);
hold on;
plot(t, all_vort(:,2), 'LineWidth', 2);
plot(t, all_vort(:,3), 'LineWidth', 2);
grid on;
ylabel('Mean vorticity');
title('Mean Vorticity');

nexttile;
plot(t, 100 * all_recirc(:,1), 'LineWidth', 2);
hold on;
plot(t, 100 * all_recirc(:,2), 'LineWidth', 2);
plot(t, 100 * all_recirc(:,3), 'LineWidth', 2);
grid on;
ylabel('Return flow (%)');
title('Upward Return-Flow Fraction');

nexttile;
plot(t, 100 * all_topCurve(:,1), 'LineWidth', 2);
hold on;
plot(t, 100 * all_topCurve(:,2), 'LineWidth', 2);
plot(t, 100 * all_topCurve(:,3), 'LineWidth', 2);
grid on;
ylabel('Particles (%)');
title('Particles Reaching Upper Inside Curve');

nexttile;
plot(t, 100 * all_exchangeIn(:,1), 'LineWidth', 2);
hold on;
plot(t, 100 * all_exchangeIn(:,2), 'LineWidth', 2);
plot(t, 100 * all_exchangeIn(:,3), 'LineWidth', 2);
grid on;
ylabel('Exchange ingress (%)');
xlabel('Time');
title('Exchange-Reservoir Particles Entering Blade Guard Region');

disp('Done.');

%% Local functions

function [fluid, solid, bellSolid, openingFluid, exchangeReservoir] = ...
    build_domain_mask_hybrid(model, X, Y, Z, zmin, zmax, wall_clearance)
    R  = sqrt(X.^2 + Y.^2);
    TH = atan2(Y, X);
    bell_zmin = -1.20;
    bell_zmax = -0.20;
    zloc = Z + 1.20;
    baseWallR = sqrt(max(1 - zloc.^2, 0));

    switch model
        case 1
            rWall = baseWallR;

        case 2
            rWall = baseWallR;

        case 3
            s = min(max((-Z - 0.20) / 1.0, 0), 1);
            rWall = baseWallR .* (1 - 0.10 * (1 + cos(6 * TH)) .* s);

        otherwise
            rWall = baseWallR;
    end

    rFluidLimit = max(rWall - wall_clearance, 0);

    %% Main internal blade guard fluid
    insideBellCurve = ...
        Z >= bell_zmin & ...
        Z <= bell_zmax & ...
        R < rFluidLimit;

    %% Fluid below blade guard opening
    lipRadius = 0.94;
    belowBell = ...
        Z >= zmin & ...
        Z < bell_zmin & ...
        R < lipRadius - wall_clearance;

    %% Exchange reservoirs
    exchangeReservoir = false(size(X));

    if model == 2
        % Model 2: side reservoir around side slits
        outerSideR = 1.14;
        exchangeReservoir = ...
            Z > -1.08 & ...
            Z < -0.66 & ...
            R > rWall + wall_clearance & ...
            R < outerSideR;

    elseif model == 3
        % Model 3: lower reservoir around daisy/lobed bottom rim
        outerSideR = 1.15;
        exchangeReservoir = ...
            Z > -1.30 & ...
            Z < -0.94 & ...
            R > 0.86 & ...
            R < outerSideR;
    end

    %% Openings through blade guard wall
    openingFluid = false(size(X));

    if model == 2
        % Four side slits
        hole_centers = [0, pi/2, pi, 3*pi/2];
        for hc = hole_centers
            dth = atan2(sin(TH - hc), cos(TH - hc));
            thisSlit = ...
                abs(dth) < 0.40 & ...
                abs(Z + 0.88) < 0.14 & ...
                abs(R - rWall) < 0.095;
            openingFluid = openingFluid | thisSlit;
        end

    elseif model == 3
        % Six exchanged openings for the daisy model
        lobe_centers = linspace(0, 2*pi, 7);
        lobe_centers(end) = [];
        for hc = lobe_centers
            dth = atan2(sin(TH - hc), cos(TH - hc));
            thisBottomOpening = ...
                abs(dth) < 0.30 & ...
                Z > -1.26 & ...
                Z < -1.02 & ...
                R > 0.72 & ...
                R < 1.08;
            openingFluid = openingFluid | thisBottomOpening;
        end
    end

    %% Fluid domain
    fluid = insideBellCurve | belowBell | exchangeReservoir | openingFluid;
    fluid = fluid & Z <= zmax;

    %% Solid blade guard wall
    shell_t = 0.055;
    bellWallCurve = ...
        Z >= bell_zmin & ...
        Z <= bell_zmax & ...
        abs(R - rWall) <= shell_t;

    if model == 2 || model == 3
        bellWallCurve(openingFluid) = false;
    end

    %% Lower side wall below blade guard opening
    lowerSideWall = ...
        Z >= zmin & ...
        Z < bell_zmin & ...
        abs(R - lipRadius) <= shell_t;

    %% Remove daisy bottom openings from lower side wall
    if model == 3
        lowerSideWall(openingFluid) = false;
    end

    %% Outer reservoir wall
    outerReservoirWall = false(size(X));

    if model == 2
        outerSideR = 1.14;
        outerReservoirWall = ...
            Z > -1.08 & ...
            Z < -0.66 & ...
            abs(R - outerSideR) <= shell_t;

    elseif model == 3
        outerSideR = 1.15;
        outerReservoirWall = ...
            Z > -1.30 & ...
            Z < -0.94 & ...
            abs(R - outerSideR) <= shell_t;
    end

    %% Central shaft
    shaft = ...
        R < 0.055 & ...
        Z > -0.48 & ...
        Z < -0.20;

    bellSolid = bellWallCurve | lowerSideWall | outerReservoirWall | shaft;
    solid = ~fluid | bellSolid;

    %% solid boundaries
    solid(1,:,:)   = true;
    solid(end,:,:) = true;
    solid(:,1,:)   = true;
    solid(:,end,:) = true;
    solid(:,:,1)   = true;
    solid(:,:,end) = true;
    fluid(solid) = false;
    openingFluid = openingFluid & fluid;
    exchangeReservoir = exchangeReservoir & fluid;
end

function [Fx, Fy, Fz] = blade_forcing_hybrid(X, Y, Z, blade_z, blade_r, omega, model)
    R  = sqrt(X.^2 + Y.^2) + 1e-8;
    TH = atan2(Y, X);

    %% Blade swirl and downward suction
    axialCore = exp(-((Z - blade_z) / 0.085).^2);
    Ftheta = omega * (R / blade_r) ...
        .* exp(-(R / (blade_r * 1.15)).^2) ...
        .* axialCore;
    Fr_blade = 13.0 ...
        .* exp(-((Z - blade_z) / 0.08).^2) ...
        .* exp(-((R - blade_r) / 0.16).^2);
    Fz_down = -15.0 ...
        .* exp(-((Z - (blade_z - 0.10)) / 0.16).^2) ...
        .* exp(-(R / (blade_r * 1.35)).^2);

    %% Internal upward wall return
    wallReturn = exp(-((R - 0.70) / 0.18).^2) ...
               .* exp(-((Z + 0.73) / 0.42).^2);
    Fz_wall_up = 14.0 * wallReturn;

    %% Upper inside-curve roll-over
    upperInsideCurve = exp(-((R - 0.62) / 0.20).^2) ...
                     .* exp(-((Z + 0.30) / 0.14).^2);
    Fr_upper_in = -8.0 * upperInsideCurve;
    Fz_upper_up = 4.2  * upperInsideCurve;

    %% Central downward return
    centralReturn = exp(-(R / 0.38).^2) ...
                  .* exp(-((Z + 0.30) / 0.22).^2);
    Fz_center_down = -5.5 * centralReturn;

    %% Lower outward spread
    lowerSpread = exp(-((Z + 1.08) / 0.18).^2) ...
                .* exp(-(R / 0.75).^2);
    Fr_lower_out = 5.0 * lowerSpread;

    %% Base cylindrical components
    Fr = Fr_blade + Fr_upper_in + Fr_lower_out;
    Fz = Fz_down + Fz_wall_up + Fz_upper_up + Fz_center_down;

    %% Model 2: side-slits
    if model == 2
        hole_centers = [0, pi/2, pi, 3*pi/2];
        slitZone = zeros(size(X));

        for hc = hole_centers
            dth = atan2(sin(TH - hc), cos(TH - hc));
            slitZone = slitZone + ...
                exp(-(dth / 0.34).^2) ...
                .* exp(-((Z + 0.88) / 0.16).^2) ...
                .* exp(-((R - 0.96) / 0.16).^2);
        end
        Fr_slit_in = -8.0 * slitZone;
        Fz_slit_mix = 3.0 * slitZone .* sin(4 * TH);
        suctionRelief = exp(-(R / 0.36).^2) ...
                      .* exp(-((Z + 0.82) / 0.20).^2);
        Fz_suction_relief = 4.0 * suctionRelief;
        Fr = Fr + Fr_slit_in;
        Fz = Fz + Fz_slit_mix + Fz_suction_relief;
    end

    %% Model 3: daisy shape with bottom exchange
    if model == 3
        lobe_centers = linspace(0, 2*pi, 7);
        lobe_centers(end) = [];
        bottomExchangeZone = zeros(size(X));

        for hc = lobe_centers
            dth = atan2(sin(TH - hc), cos(TH - hc));
            bottomExchangeZone = bottomExchangeZone + ...
                exp(-(dth / 0.28).^2) ...
                .* exp(-((Z + 1.12) / 0.16).^2) ...
                .* exp(-((R - 0.96) / 0.18).^2);
        end

        Fr_bottom_in = -6.5 * bottomExchangeZone;
        Fz_bottom_up = 5.5 * bottomExchangeZone;
        Fz_lobed_mix = 2.4 * bottomExchangeZone .* sin(6 * TH);

        daisySuctionRelief = exp(-(R / 0.42).^2) ...
                           .* exp(-((Z + 1.02) / 0.22).^2);

        Fz_daisy_relief = 3.2 * daisySuctionRelief;
        Fr = Fr + Fr_bottom_in;
        Fz = Fz + Fz_bottom_up + Fz_lobed_mix + Fz_daisy_relief;
    end

    Fx = Fr .* cos(TH) - Ftheta .* sin(TH);
    Fy = Fr .* sin(TH) + Ftheta .* cos(TH);
end

function phi_new = advect_field(phi, U, V, W, X, Y, Z, dt, fluid)
    Xb = X - dt * U;
    Yb = Y - dt * V;
    Zb = Z - dt * W;
    phi_new = interpn(X, Y, Z, phi, Xb, Yb, Zb, 'linear', 0);
    phi_new(~fluid) = 0;
end

function L = lap3d(A, dx, dy, dz)
    L = ...
        (circshift(A, [-1  0  0]) - 2*A + circshift(A, [1 0 0])) / dx^2 + ...
        (circshift(A, [ 0 -1  0]) - 2*A + circshift(A, [0 1 0])) / dy^2 + ...
        (circshift(A, [ 0  0 -1]) - 2*A + circshift(A, [0 0 1])) / dz^2;
end

function [U, V, W] = apply_noslip(U, V, W, solid)
    U(solid) = 0;
    V(solid) = 0;
    W(solid) = 0;
    U(1,:,:) = 0;   U(end,:,:) = 0;
    U(:,1,:) = 0;   U(:,end,:) = 0;
    U(:,:,1) = 0;   U(:,:,end) = 0;
    V(1,:,:) = 0;   V(end,:,:) = 0;
    V(:,1,:) = 0;   V(:,end,:) = 0;
    V(:,:,1) = 0;   V(:,:,end) = 0;
    W(1,:,:) = 0;   W(end,:,:) = 0;
    W(:,1,:) = 0;   W(:,end,:) = 0;
    W(:,:,1) = 0;   W(:,:,end) = 0;
end

function [U, V, W, P] = pressure_project(U, V, W, P, fluid, h, dt, nIter)
    div = divergence3d(U, V, W, h, h, h);
    div(~fluid) = 0;

    for it = 1:nIter

        Pe = circshift(P, [-1  0  0]);
        Pw = circshift(P, [ 1  0  0]);
        Pn = circshift(P, [ 0 -1  0]);
        Ps = circshift(P, [ 0  1  0]);
        Pt = circshift(P, [ 0  0 -1]);
        Pb = circshift(P, [ 0  0  1]);

        Fe = circshift(fluid, [-1  0  0]);
        Fw = circshift(fluid, [ 1  0  0]);
        Fn = circshift(fluid, [ 0 -1  0]);
        Fs = circshift(fluid, [ 0  1  0]);
        Ft = circshift(fluid, [ 0  0 -1]);
        Fb = circshift(fluid, [ 0  0  1]);

        sumP = Pe.*Fe + Pw.*Fw + Pn.*Fn + Ps.*Fs + Pt.*Ft + Pb.*Fb;
        cnt  = double(Fe) + double(Fw) + double(Fn) + double(Fs) + double(Ft) + double(Fb);
        cnt(cnt == 0) = 1;
        Pnew = (sumP - (h^2) * (div / dt)) ./ cnt;
        Pnew(~fluid) = 0;
        P = Pnew;
    end

    [dPdx, dPdy, dPdz] = grad3d(P, h, h, h);
    U = U - dt * dPdx;
    V = V - dt * dPdy;
    W = W - dt * dPdz;
    U(~fluid) = 0;
    V(~fluid) = 0;
    W(~fluid) = 0;
end

function div = divergence3d(U, V, W, dx, dy, dz)
    dUdx = (circshift(U, [-1 0 0]) - circshift(U, [1 0 0])) / (2 * dx);
    dVdy = (circshift(V, [0 -1 0]) - circshift(V, [0 1 0])) / (2 * dy);
    dWdz = (circshift(W, [0 0 -1]) - circshift(W, [0 0 1])) / (2 * dz);
    div = dUdx + dVdy + dWdz;
end

function [dAdx, dAdy, dAdz] = grad3d(A, dx, dy, dz)
    dAdx = (circshift(A, [-1 0 0]) - circshift(A, [1 0 0])) / (2 * dx);
    dAdy = (circshift(A, [0 -1 0]) - circshift(A, [0 1 0])) / (2 * dy);
    dAdz = (circshift(A, [0 0 -1]) - circshift(A, [0 0 1])) / (2 * dz);
end

function [omx, omy, omz, ommag] = curl3d(U, V, W, dx, dy, dz)
    [~,    dUdy, dUdz] = grad3d(U, dx, dy, dz);
    [dVdx, ~,    dVdz] = grad3d(V, dx, dy, dz);
    [dWdx, dWdy, ~   ] = grad3d(W, dx, dy, dz);
    omx = dWdy - dVdz;
    omy = dUdz - dWdx;
    omz = dVdx - dUdy;
    ommag = sqrt(omx.^2 + omy.^2 + omz.^2);
end

function [px, py, pz, ptag] = seed_particles_hybrid(Np, fluid, exchangeReservoir, X, Y, Z, model)
    R = sqrt(X.^2 + Y.^2);
    lowerFood = fluid & ...
        Z < -1.05 & ...
        R < 0.85;
    bladeRegion = fluid & ...
        Z > -0.75 & Z < -0.35 & ...
        R < 0.75;
    upperInsideBell = fluid & ...
        Z > -0.45 & Z < -0.22 & ...
        R > 0.30 & R < 0.88;
    exchangeFood = exchangeReservoir;

    if model == 2
        nExchange = round(0.25 * Np);
        nLower    = round(0.35 * Np);
        nBlade    = round(0.25 * Np);
        nUpper    = Np - nExchange - nLower - nBlade;

    elseif model == 3
        nExchange = round(0.25 * Np);
        nLower    = round(0.35 * Np);
        nBlade    = round(0.25 * Np);
        nUpper    = Np - nExchange - nLower - nBlade;

    else
        nExchange = 0;
        nLower    = round(0.45 * Np);
        nBlade    = round(0.35 * Np);
        nUpper    = Np - nLower - nBlade;
    end

    idxLower    = find(lowerFood);
    idxBlade    = find(bladeRegion);
    idxUpper    = find(upperInsideBell);
    idxExchange = find(exchangeFood);

    if isempty(idxLower),    idxLower    = find(fluid); end
    if isempty(idxBlade),    idxBlade    = find(fluid); end
    if isempty(idxUpper),    idxUpper    = find(fluid); end
    if isempty(idxExchange), idxExchange = find(fluid); end

    pickLower = idxLower(randi(numel(idxLower), [nLower 1]));
    pickBlade = idxBlade(randi(numel(idxBlade), [nBlade 1]));
    pickUpper = idxUpper(randi(numel(idxUpper), [nUpper 1]));
    pick = [pickLower; pickBlade; pickUpper];
    ptag = ones(numel(pick), 1);

    if nExchange > 0
        pickExchange = idxExchange(randi(numel(idxExchange), [nExchange 1]));
        pick = [pick; pickExchange];
        ptag = [ptag; 2 * ones(nExchange, 1)];
    end

    px = X(pick);
    py = Y(pick);
    pz = Z(pick);
end

function [px, py, pz] = advect_particles(px, py, pz, U, V, W, X, Y, Z, dt, fluid)

    Ui = interpn(X, Y, Z, U, px, py, pz, 'linear', 0);
    Vi = interpn(X, Y, Z, V, px, py, pz, 'linear', 0);
    Wi = interpn(X, Y, Z, W, px, py, pz, 'linear', 0);

    %% Unresolved turbulent food mixing
    turb = 0.018;

    jx = turb * randn(size(px));
    jy = turb * randn(size(py));
    jz = 0.65 * turb * randn(size(pz));
    pxt = px + dt * Ui + jx * sqrt(dt);
    pyt = py + dt * Vi + jy * sqrt(dt);
    pzt = pz + dt * Wi + jz * sqrt(dt);

    valid = interpn(X, Y, Z, double(fluid), pxt, pyt, pzt, 'nearest', 0) > 0.5;

    px(valid) = pxt(valid);
    py(valid) = pyt(valid);
    pz(valid) = pzt(valid);

    %% Wall rebound approximation
    hit = ~valid;

    if any(hit)
        px(hit) = px(hit) - 0.40 * dt * Ui(hit) + 0.006 * randn(sum(hit), 1);
        py(hit) = py(hit) - 0.40 * dt * Vi(hit) + 0.006 * randn(sum(hit), 1);
        pz(hit) = pz(hit) - 0.40 * dt * Wi(hit) + 0.004 * randn(sum(hit), 1);
    end

    stillValid = interpn(X, Y, Z, double(fluid), px, py, pz, 'nearest', 0) > 0.5;

    if any(~stillValid)
        idxFluid = find(fluid);
        bad = find(~stillValid);
        replacement = idxFluid(randi(numel(idxFluid), [numel(bad) 1]));
        px(bad) = X(replacement);
        py(bad) = Y(replacement);
        pz(bad) = Z(replacement);
    end
end

function [Xs, Ys, Zs] = bell_surface(model, res)
    [U, V] = meshgrid(linspace(0, 2*pi, res), linspace(0, pi/2, res));

    switch model

        case 1
            Xs = sin(V) .* cos(U);
            Ys = sin(V) .* sin(U);
            Zs = -(1.2 - cos(V));
            mask = true(size(Xs));

        case 2
            Xs = sin(V) .* cos(U);
            Ys = sin(V) .* sin(U);
            Zs = -(1.2 - cos(V));
            mask = true(size(Xs));
            hole_centers = [0, pi/2, pi, 3*pi/2];

            for hc = hole_centers
                dth = atan2(sin(U - hc), cos(U - hc));
                mask((abs(dth) < 0.45) & (abs(V - 1.25) < 0.12)) = false;
            end

        case 3
            Rmod = 1.0 - 0.12 * (1 + cos(6 * U)) .* (V / (pi/2));
            Xs = Rmod .* sin(V) .* cos(U);
            Ys = Rmod .* sin(V) .* sin(U);
            Zs = -(1.2 - cos(V));
            mask = V < (pi/2 - 0.25 * (1 + cos(6 * U)) / 2);
    end

    Xs(~mask) = nan;
    Ys(~mask) = nan;
    Zs(~mask) = nan;
end

function [xline, zline] = bell_outline_xz(model, n)
    v = linspace(0, pi/2, n);
    zline = -(1.2 - cos(v));
    r = sin(v);

    switch model

        case 1
            xline = r;

        case 2
            xline = r;

        case 3
            u0 = 0;
            Rmod = 1.0 - 0.12 * (1 + cos(6 * u0)) .* (v / (pi/2));
            xline = Rmod .* r;
    end
end

function draw_inside_top_curve_region()
    th = linspace(0, 2*pi, 140);
    z1 = -0.45;
    z2 = -0.22;
    r1 = 0.45;
    r2 = 0.88;
    [T, ZZ] = meshgrid(th, linspace(z1, z2, 8));
    X1 = r1 * cos(T);
    Y1 = r1 * sin(T);
    X2 = r2 * cos(T);
    Y2 = r2 * sin(T);

    surf(X1, Y1, ZZ, ...
        'FaceColor', [1.0 0.8 0.1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.055);

    surf(X2, Y2, ZZ, ...
        'FaceColor', [1.0 0.8 0.1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.055);
end

function draw_exchange_reservoir_region(model)
    th = linspace(0, 2*pi, 140);

    if model == 2
        z1 = -1.08;
        z2 = -0.66;
        r1 = 0.98;
        r2 = 1.14;

    elseif model == 3
        z1 = -1.30;
        z2 = -0.94;
        r1 = 0.86;
        r2 = 1.15;

    else
        return;
    end

    [T, ZZ] = meshgrid(th, linspace(z1, z2, 8));
    X1 = r1 * cos(T);
    Y1 = r1 * sin(T);
    X2 = r2 * cos(T);
    Y2 = r2 * sin(T);

    surf(X1, Y1, ZZ, ...
        'FaceColor', [0.1 0.9 1.0], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.045);

    surf(X2, Y2, ZZ, ...
        'FaceColor', [0.1 0.9 1.0], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.045);
end

function m = safe_mean(A)
    A = A(:);

    if isempty(A)
        m = 0;
    else
        m = mean(A);
    end
end

%% Design evaluation

fprintf('Design evaluation\n');
startIndex = round(0.40 * nSteps);

avgSpeed    = mean(all_speed(startIndex:end, :), 1);
avgVort     = mean(all_vort(startIndex:end, :), 1);
avgRecirc   = mean(all_recirc(startIndex:end, :), 1);
avgExchange = mean(all_exchangeIn(startIndex:end, :), 1);
avgTopCurve = mean(all_topCurve(startIndex:end, :), 1);

normSpeed    = normalize_metric(avgSpeed);
normVort     = normalize_metric(avgVort);
normRecirc   = normalize_metric(avgRecirc);
normExchange = normalize_metric(avgExchange);
normTopCurve = normalize_metric(avgTopCurve);

suctionLockRisk = ...
    (0.45 * normSpeed + ...
     0.45 * normVort + ...
     0.10 * normRecirc) ...
     .* (1 - normExchange);

antiSuctionScore = 1 - suctionLockRisk;

usefulMixingScore = ...
    0.18 * normSpeed + ...
    0.18 * normVort + ...
    0.18 * normRecirc + ...
    0.36 * normExchange + ...
    0.10 * normTopCurve;

efficiencyScore = ...
    0.70 * usefulMixingScore + ...
    0.30 * antiSuctionScore;

%%   real mixing = exchange ingress × vorticity × return flow

realMixingIndex = avgExchange .* avgVort .* avgRecirc;
if max(realMixingIndex) > 0
    normalizedRealMixingIndex = realMixingIndex ./ max(realMixingIndex);
else
    normalizedRealMixingIndex = zeros(size(realMixingIndex));
end

fprintf('Detailed values\n');

for i = 1:nModels
    fprintf('%s\n', model_names{i});
    fprintf('  Mean speed:                 %.4f\n', avgSpeed(i));
    fprintf('  Mean vorticity:             %.4f\n', avgVort(i));
    fprintf('  Return flow:                %.2f %%\n', 100 * avgRecirc(i));
    fprintf('  Exchange ingress:           %.2f %%\n', 100 * avgExchange(i));
    fprintf('  Upper-curve reach:          %.2f %%\n', 100 * avgTopCurve(i));
    fprintf('  Suction-lock risk:          %.3f / 1.000\n', suctionLockRisk(i));
    fprintf('  Anti-suction score:         %.3f / 1.000\n', antiSuctionScore(i));
    fprintf('  Useful mixing score:        %.3f / 1.000\n', usefulMixingScore(i));
    fprintf('  Final efficiency score:     %.3f / 1.000\n', efficiencyScore(i));
    fprintf('  Real mixing index:          %.4f\n', realMixingIndex(i));
    fprintf('  Normalized real mixing:     %.3f / 1.000\n\n', normalizedRealMixingIndex(i));
end

[bestEfficiencyScore, bestEfficiencyIndex] = max(efficiencyScore);
[bestRealMixingScore, bestRealMixingIndex] = max(realMixingIndex);

fprintf('Evaluation Ranking\n');

fprintf('Best model by improved weighted efficiency score:\n');
fprintf('  %s\n', model_names{bestEfficiencyIndex});
fprintf('  Score: %.3f / 1.000\n\n', bestEfficiencyScore);

fprintf('Best model by real mixing index:\n');
fprintf('  %s\n', model_names{bestRealMixingIndex});
fprintf('  Real mixing index: %.4f\n', bestRealMixingScore);
fprintf('  Normalized score: %.3f / 1.000\n\n', normalizedRealMixingIndex(bestRealMixingIndex));

%% printing full rankings

[~, rankEfficiency] = sort(efficiencyScore, 'descend');
[~, rankRealMixing] = sort(realMixingIndex, 'descend');

fprintf('Ranking by improved weighted efficiency score:\n');
for k = 1:nModels
    idx = rankEfficiency(k);
    fprintf('  %d. %s  —  %.3f / 1.000\n', ...
        k, model_names{idx}, efficiencyScore(idx));
end

fprintf('\nRanking by real mixing index:\n');
for k = 1:nModels
    idx = rankRealMixing(k);
    fprintf('  %d. %s  —  %.4f\n', ...
        k, model_names{idx}, realMixingIndex(idx));
end

%% plotting improved efficiency score

figure('Color','w','Name','Improved Efficiency Score','Position',[260 180 1000 520]);

bar(efficiencyScore);
set(gca, ...
    'XTick', 1:nModels, ...
    'XTickLabel', model_names, ...
    'XTickLabelRotation', 20);

ylabel('Improved weighted efficiency score');
title('Overall Blade Guard Efficiency Score with Suction-Lock Penalty');
ylim([0, max(1.0, 1.10 * max(efficiencyScore))]);
grid on;

for i = 1:nModels
    text(i, efficiencyScore(i) + 0.03, ...
        sprintf('%.3f', efficiencyScore(i)), ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold');
end

%% plotting real mixing index

figure('Color','w','Name','Real Mixing Index','Position',[300 220 1000 520]);

bar(normalizedRealMixingIndex);
set(gca, ...
    'XTick', 1:nModels, ...
    'XTickLabel', model_names, ...
    'XTickLabelRotation', 20);

ylabel('Normalized real mixing index');
title('Real Mixing Index: Exchange × Vorticity × Return Flow');
ylim([0, max(1.0, 1.10 * max(normalizedRealMixingIndex))]);
grid on;

for i = 1:nModels
    text(i, normalizedRealMixingIndex(i) + 0.03, ...
        sprintf('%.3f', normalizedRealMixingIndex(i)), ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold');
end

%% plotting evaluation metrics

figure('Color','w','Name','Normalized Design Metrics','Position',[240 160 1200 620]);

metricMatrix = [
    normSpeed(:), ...
    normVort(:), ...
    normRecirc(:), ...
    normExchange(:), ...
    normTopCurve(:), ...
    antiSuctionScore(:), ...
    normalizedRealMixingIndex(:)
];

bar(metricMatrix);

set(gca, ...
    'XTick', 1:nModels, ...
    'XTickLabel', model_names, ...
    'XTickLabelRotation', 20);

ylabel('Normalized value');
title('Normalized Metrics Used to Compare Guard Designs');
legend({ ...
    'Speed', ...
    'Vorticity', ...
    'Return flow', ...
    'Exchange ingress', ...
    'Upper curve reach', ...
    'Anti-suction', ...
    'Real mixing index'}, ...
    'Location','eastoutside');

grid on;

%% metric evaluation normalization

function y = normalize_metric(x)
    xmin = min(x);
    xmax = max(x);

    if abs(xmax - xmin) < 1e-12
        y = ones(size(x));
    else
        y = (x - xmin) ./ (xmax - xmin);
    end
end