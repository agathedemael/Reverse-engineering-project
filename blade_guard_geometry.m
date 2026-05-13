% Blade Guard Geometry Only
clear; clc; close all;

res = 60;
[U_grid, V_grid] = meshgrid(linspace(0, 2*pi, res), linspace(0, pi/2, res));
[U_o, V_o] = meshgrid(linspace(0, 2*pi, res), linspace(0, 1.4, res));

model_names = {'1. Base model (closed)', '2. Perforated (wide slits)', ...
               '3. Daisy shape (slightly recessed)'};

figure('Color', 'w', 'Name', 'Blade guard geometry');

for i = 1:3
    subplot(2, 2, i); hold on;

    switch i
        case 1 % Baseline
            X = sin(V_grid).*cos(U_grid);
            Y = sin(V_grid).*sin(U_grid);
            Z = 1.2 - cos(V_grid);
            mask = true(size(X));

        case 2 % Perforated (reduced height holes)
            X = sin(V_grid).*cos(U_grid);
            Y = sin(V_grid).*sin(U_grid);
            Z = 1.2 - cos(V_grid);
            mask = true(size(X));

            hole_centers = [0, pi/2, pi, 3*pi/2];
            for hc = hole_centers
                mask((abs(U_grid-hc) < 0.45) & (abs(V_grid-1.25) < 0.10)) = false;
            end

        case 3 % Subtle Daisy (deeper recesses)
            R_mod = 1.0 - 0.12 * (1 + cos(6*U_grid)) .* (V_grid/(pi/2));
            X = R_mod .* sin(V_grid).*cos(U_grid);
            Y = R_mod .* sin(V_grid).*sin(U_grid);
            Z = 1.2 - cos(V_grid);

            % Cutoff mask also using 6 peaks
            mask = V_grid < (pi/2 - 0.25 * (1 + cos(6*U_grid))/2);
    end

    Xp = X; Yp = Y; Zp = Z;
    Xp(~mask) = nan;

    surf(Xp, Yp, Zp, ...
        'FaceColor', [0.4 0.7 0.9], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.25);

    camlight; lighting phong;
    title(model_names{i});
    view(35, 25); axis equal; axis off;
end

rotate3d on;

function [Xp, Yp, Zp] = bell_geometry(i, res)
    [U_grid, V_grid] = meshgrid(linspace(0, 2*pi, res), linspace(0, pi/2, res));

    switch i
        case 1 % Baseline
            X = sin(V_grid).*cos(U_grid);
            Y = sin(V_grid).*sin(U_grid);
            Z = 1.2 - cos(V_grid);
            mask = true(size(X));

        case 2 % Perforated
            X = sin(V_grid).*cos(U_grid);
            Y = sin(V_grid).*sin(U_grid);
            Z = 1.2 - cos(V_grid);
            mask = true(size(X));

            hole_centers = [0, pi/2, pi, 3*pi/2];
            for hc = hole_centers
                mask((abs(U_grid-hc) < 0.45) & (abs(V_grid-1.25) < 0.10)) = false;
            end

        case 3 % Subtle Daisy
            R_mod = 1.0 - 0.12 * (1 + cos(6*U_grid)) .* (V_grid/(pi/2));
            X = R_mod .* sin(V_grid).*cos(U_grid);
            Y = R_mod .* sin(V_grid).*sin(U_grid);
            Z = -(1.2 - cos(V_grid));

            mask = V_grid < (pi/2 - 0.25 * (1 + cos(6*U_grid))/2);
    end

    Xp = X; Yp = Y; Zp = Z;
    Xp(~mask) = nan;
    Yp(~mask) = nan;
    Zp(~mask) = nan;
end

function write_stl_ascii(filename, TR)
    F = TR.ConnectivityList;
    V = TR.Points;

    fid = fopen(filename, 'w');
    fprintf(fid, 'solid bladeguard\n');

    for i = 1:size(F,1)
        v1 = V(F(i,1),:);
        v2 = V(F(i,2),:);
        v3 = V(F(i,3),:);

        n = cross(v2 - v1, v3 - v1);
        if norm(n) > 0
            n = n / norm(n);
        else
            n = [0 0 0];
        end

        fprintf(fid, ' facet normal %.6f %.6f %.6f\n', n);
        fprintf(fid, '  outer loop\n');
        fprintf(fid, '   vertex %.6f %.6f %.6f\n', v1);
        fprintf(fid, '   vertex %.6f %.6f %.6f\n', v2);
        fprintf(fid, '   vertex %.6f %.6f %.6f\n', v3);
        fprintf(fid, '  endloop\n');
        fprintf(fid, ' endfacet\n');
    end

    fprintf(fid, 'endsolid bladeguard\n');
    fclose(fid);
end

function export_bell_stl()
    res = 80;
    model_names = {'baseline', 'perforated', 'daisy'};

    for i = 1:3
        [X, Y, Z] = bell_geometry(i, res);

        % --- 1. Build UV grid ---
        [U, V] = meshgrid(linspace(0,1,size(X,2)), linspace(0,1,size(X,1)));

        % --- 2. Remove NaNs ---
        valid = ~isnan(X);
        Uv = U(valid);
        Vv = V(valid);
        Xv = X(valid);
        Yv = Y(valid);
        Zv = Z(valid);

        % --- 3. Remove duplicate UV points ---
        UV = [Uv(:), Vv(:)];
        [UV_unique, ia, ~] = unique(UV, 'rows');
        Xv = Xv(ia);
        Yv = Yv(ia);
        Zv = Zv(ia);

        % --- 4. Delaunay triangulation in UV ---
        DT = delaunayTriangulation(UV_unique);

        % --- 5. Build triangulation in XYZ ---
        TR = triangulation(DT.ConnectivityList, [Xv, Yv, Zv]);

        % --- 6. Export STL using ASCII writer ---
        write_stl_ascii([model_names{i} '.stl'], TR);

        fprintf("Exported %s.stl with %d vertices and %d faces\n", ...
            model_names{i}, size(TR.Points,1), size(TR.ConnectivityList,1));
    end
end

export_bell_stl();

