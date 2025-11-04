%% analyzeWaterRetainingWallCSA
% MATLAB script to evaluate overturning, sliding, and bearing capacity for a
% cantilever water-retaining wall using CSA-style limit state checks.
%
% All input dimensions are in metres and loads are per metre run of wall.
% The script is organised so that geometry, materials, and load assumptions
% are declared once and can be quickly updated to match project-specific
% details.  Where the specification was silent (e.g., stem thickness),
% typical starting values are provided and clearly flagged for update.
%
% Outputs:
%   * Factors of safety against overturning and sliding
%   * Bearing pressure distribution and comparison with factored soil
%     resistance
%   * Load breakdown tables for traceability

clear; clc;

input = getDefaultInput();

loads = computeLoads(input);
results = evaluateStability(loads, input);

reportLoads(loads);
reportResults(results);

%% --- Helper Functions --------------------------------------------------

function input = getDefaultInput()
    % Geometric inputs (metres)
    mm = 1e-3;

    geom = struct();
    geom.stem_height = 3305 * mm;              % From footing top to crest
    geom.stem_thickness_base = 0.30;           % Assumed base stem thickness (update as required)
    geom.stem_thickness_top  = 0.20;           % Assumed stem thickness at top (update as required)
    geom.base_thickness      = 0.45;           % Assumed footing thickness (update as required)
    geom.toe_width           = 610 * mm;
    geom.heel_width          = 1270 * mm;
    geom.shear_key_width     = 280 * mm;
    geom.shear_key_depth     = 370 * mm;       % Extent below footing soffit
    geom.shear_key_offset    = 0.0;            % Offset from toe edge (0 => flush with toe)

    % Gutter geometry (relative to stem)
    gutter = struct();
    gutter.depth_from_top        = 400 * mm;   % Vertical offset from stem top to gutter soffit
    gutter.horizontal_projection = 800 * mm;   % Projection from stem towards heel
    gutter.horizontal_thickness  = 150 * mm;
    gutter.vertical_height       = 290 * mm;
    gutter.vertical_thickness    = 150 * mm;
    gutter.orientation           = "heel";     % Use "toe" if gutter projects towards the toe

    % Roof slab (tributary width supplied)
    roof = struct();
    roof.thickness = 225 * mm;
    roof.tributary_span = 2.75;               % Effective width contributing to wall

    % Materials
    materials = struct();
    materials.gamma_concrete = 24.0;          % kN/m^3
    materials.gamma_water    = 9.81;          % kN/m^3
    materials.gamma_soil     = 19.0;          % kN/m^3 (foundation soil)
    materials.phi_soil       = deg2rad(32);   % Friction angle of supporting soil
    materials.delta_interface = deg2rad(15);  % Concrete-soil interface friction angle

    % Loading
    loading = struct();
    loading.water_depth          = 1.60;      % metres above heel
    loading.include_heel_water_weight = true; % Adds vertical water load over heel slab
    loading.uplift_reduction_factor = 1.0;    % 1.0 = full hydrostatic uplift under heel

    % Design/limit state factors (CSA-inspired)
    design = struct();
    design.FS_overturning = 1.50;             % Minimum resisting/demand ratio
    design.FS_sliding     = 1.50;
    design.resistance_factor_bearing = 0.50;  % φ_g for factored bearing resistance
    design.foundation_depth = 0.00;           % Embedment depth of footing underside (m)

    input = struct();
    input.geometry = geom;
    input.geometry.gutter = gutter;
    input.geometry.roof   = roof;
    input.materials = materials;
    input.loading    = loading;
    input.design     = design;
end

function loads = computeLoads(input)
    g  = input.geometry;
    m  = input.materials;
    ld = input.loading;

    vertical = struct('name',{},'value',{},'sign',{},'x',{},'z',{});

    % --- Footing self-weight broken into toe, stem block, and heel segments
    vertical = addVerticalLoad(vertical, "Toe footing self-weight", ...
        m.gamma_concrete * g.toe_width * g.base_thickness, "down", g.toe_width/2, g.base_thickness/2);

    vertical = addVerticalLoad(vertical, "Stem block footing weight", ...
        m.gamma_concrete * g.stem_thickness_base * g.base_thickness, "down", ...
        g.toe_width + g.stem_thickness_base/2, g.base_thickness/2);

    vertical = addVerticalLoad(vertical, "Heel footing self-weight", ...
        m.gamma_concrete * g.heel_width * g.base_thickness, "down", ...
        g.toe_width + g.stem_thickness_base + g.heel_width/2, g.base_thickness/2);

    % --- Stem self-weight (tapered section)
    stem_area = 0.5 * (g.stem_thickness_base + g.stem_thickness_top) * g.stem_height;
    stem_weight = m.gamma_concrete * stem_area;
    stem_centroid_z = g.base_thickness + centroidZTrapezoid(g.stem_height, ...
        g.stem_thickness_base, g.stem_thickness_top);

    vertical = addVerticalLoad(vertical, "Stem concrete weight", stem_weight, ...
        "down", g.toe_width + g.stem_thickness_base/2, stem_centroid_z);

    % --- Roof slab weight (tributary loading)
    roof = g.roof;
    roof_weight = m.gamma_concrete * roof.thickness * roof.tributary_span;
    roof_centroid_z = g.base_thickness + g.stem_height + roof.thickness/2;
    vertical = addVerticalLoad(vertical, "Roof slab weight", roof_weight, ...
        "down", g.toe_width + g.stem_thickness_base/2, roof_centroid_z);

    % --- Shear key weight
    shear_key_volume = g.shear_key_width * g.shear_key_depth;
    shear_key_weight = m.gamma_concrete * shear_key_volume;
    shear_key_x = g.shear_key_offset + g.shear_key_width/2;
    shear_key_z = -g.shear_key_depth/2;
    vertical = addVerticalLoad(vertical, "Shear key weight", shear_key_weight, ...
        "down", shear_key_x, shear_key_z);

    % --- Gutter components
    gutter = g.gutter;
    if gutter.orientation == "heel"
        gutter_x_base = g.toe_width + g.stem_thickness_base;
    else
        gutter_x_base = 0.0;
    end

    % Horizontal portion
    gutter_horizontal_volume = gutter.horizontal_projection * gutter.horizontal_thickness;
    gutter_horizontal_weight = m.gamma_concrete * gutter_horizontal_volume;
    gutter_top_elevation = g.base_thickness + g.stem_height - gutter.depth_from_top;
    gutter_horizontal_z = gutter_top_elevation - gutter.horizontal_thickness/2;
    gutter_horizontal_x = gutter_x_base + gutter.horizontal_projection/2;

    vertical = addVerticalLoad(vertical, "Gutter horizontal plate weight", ...
        gutter_horizontal_weight, "down", gutter_horizontal_x, gutter_horizontal_z);

    % Vertical upstand
    gutter_vertical_volume = gutter.vertical_height * gutter.vertical_thickness;
    gutter_vertical_weight = m.gamma_concrete * gutter_vertical_volume;
    gutter_vertical_z = gutter_top_elevation + gutter.vertical_height/2;
    gutter_vertical_x = gutter_x_base + gutter.horizontal_projection - gutter.vertical_thickness/2;

    vertical = addVerticalLoad(vertical, "Gutter vertical upstand weight", ...
        gutter_vertical_weight, "down", gutter_vertical_x, gutter_vertical_z);

    % --- Retained water weight over heel
    if ld.include_heel_water_weight
        heel_water_volume = g.heel_width * ld.water_depth;
        heel_water_weight = m.gamma_water * heel_water_volume;
        heel_water_z = g.base_thickness + ld.water_depth/2;
        heel_water_x = g.toe_width + g.stem_thickness_base + g.heel_width/2;

        vertical = addVerticalLoad(vertical, "Water weight over heel", ...
            heel_water_weight, "down", heel_water_x, heel_water_z);
    end

    % --- Hydrostatic uplift beneath heel (reduction factor allows tailoring)
    if ld.uplift_reduction_factor > 0
        uplift_pressure = ld.uplift_reduction_factor * m.gamma_water * ld.water_depth;
        uplift_force = uplift_pressure * g.heel_width;
        uplift_x = g.toe_width + g.stem_thickness_base + g.heel_width/2;
        uplift_z = -g.base_thickness/2; % act at footing soffit centroid under heel

        vertical = addVerticalLoad(vertical, "Hydrostatic uplift under heel", ...
            uplift_force, "up", uplift_x, uplift_z);
    end

    % --- Hydrostatic load on stem
    hydro = struct();
    hydro.name = "Hydrostatic pressure on stem";
    hydro.force = 0.5 * m.gamma_water * ld.water_depth^2;
    hydro.resultant_z = g.base_thickness + ld.water_depth/3; % from footing bottom
    hydro.point_of_action_x = g.toe_width + g.stem_thickness_base; % heel face of stem

    loads = struct();
    loads.vertical = vertical;
    loads.hydrostatic = hydro;
end

function results = evaluateStability(loads, input)
    design = input.design;
    geom  = input.geometry;
    soil  = input.materials;

    vertical = loads.vertical;
    hydro = loads.hydrostatic;

    values = [vertical.value];
    signs  = [vertical.sign];
    arms   = [vertical.x];

    Rz = sum(values .* signs);
    Mr = sum(values .* signs .* arms);
    Mo = hydro.force * hydro.resultant_z;

    FS_overturn = Mr / Mo;

    % Sliding resistance
    H_total = hydro.force; % Driving force towards toe

    friction_resistance = max(Rz, 0) * tan(soil.delta_interface);
    passive_coefficient = tan(pi/4 + soil.phi_soil/2)^2;
    passive_resistance = 0.5 * soil.gamma_soil * geom.shear_key_depth^2 * passive_coefficient;
    sliding_resistance = friction_resistance + passive_resistance;
    FS_sliding = sliding_resistance / H_total;

    % Bearing pressures
    base_width = geom.toe_width + geom.stem_thickness_base + geom.heel_width;
    area = base_width * 1.0; % per metre length

    net_moment_toe = Mr - Mo;
    resultant_x = net_moment_toe / Rz; % distance from toe

    eccentricity = base_width/2 - resultant_x;
    q_avg = Rz / area;
    q_max = q_avg * (1 + 6*eccentricity/base_width);
    q_min = q_avg * (1 - 6*eccentricity/base_width);

    % Bearing capacity (Terzaghi general shear, cohesionless)
    phi = soil.phi_soil;
    tan_phi = tan(phi);
    Nq = exp(pi * tan_phi) * tan(pi/4 + phi/2)^2;
    Ngamma = 2 * (Nq + 1) * tan_phi;

    B_eff = max(base_width - 2*abs(eccentricity), 0.01);
    Df = design.foundation_depth;

    q_ult = 0.5 * soil.gamma_soil * B_eff * Ngamma + soil.gamma_soil * Df * Nq;
    q_factored = design.resistance_factor_bearing * q_ult;
    FS_bearing = q_factored / q_max;

    results = struct();
    results.overturning = struct('Mr',Mr,'Mo',Mo,'FS',FS_overturn,'FS_required',design.FS_overturning);
    results.sliding = struct('Rz',Rz,'friction',friction_resistance,'passive',passive_resistance, ...
        'resistance',sliding_resistance,'H',H_total,'FS',FS_sliding,'FS_required',design.FS_sliding);
    results.bearing = struct('Rz',Rz,'resultant_x',resultant_x,'eccentricity',eccentricity, ...
        'q_avg',q_avg,'q_max',q_max,'q_min',q_min,'q_factored',q_factored,'FS',FS_bearing);
end

function reportLoads(loads)
    vertical = loads.vertical;

    fprintf('\nVertical load summary (positive = downward)\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('%-40s %10s %10s\n','Component','Load (kN/m)','Lever arm (m)');
    fprintf('-------------------------------------------------------------\n');

    for k = 1:numel(vertical)
        load_sign = vertical(k).sign;
        load_value = vertical(k).value * load_sign;
        fprintf('%-40s %10.2f %10.3f\n', vertical(k).name, load_value, vertical(k).x);
    end

    hydro = loads.hydrostatic;
    fprintf('\nHydrostatic load: %.2f kN/m at z = %.3f m above footing bottom\n', ...
        hydro.force, hydro.resultant_z);
end

function reportResults(results)
    fprintf('\nCSA-style stability check summary\n');
    fprintf('=============================================================\n');

    o = results.overturning;
    fprintf('Overturning: Mr = %.2f kN-m/m, Mo = %.2f kN-m/m, FS = %.2f (>= %.2f?)\n', ...
        o.Mr, o.Mo, o.FS, o.FS_required);

    s = results.sliding;
    fprintf('Sliding: Rz = %.2f kN/m, R_fric = %.2f, R_passive = %.2f, H = %.2f, FS = %.2f (>= %.2f?)\n', ...
        s.Rz, s.friction, s.passive, s.H, s.FS, s.FS_required);

    b = results.bearing;
    fprintf('Bearing: q_avg = %.2f kPa, q_max = %.2f kPa, q_min = %.2f kPa, phi_g*q_ult = %.2f kPa, FS = %.2f\n', ...
        kPa(b.q_avg), kPa(b.q_max), kPa(b.q_min), kPa(b.q_factored), b.FS);
    fprintf('Resultant location: x = %.3f m from toe, eccentricity = %.3f m\n', ...
        b.resultant_x, b.eccentricity);

    fprintf('=============================================================\n');
    fprintf('NOTE: Review assumed stem/base thicknesses, uplift factor, and embedment depth\n');
    fprintf('      to match project-specific design per CSA requirements.\n');
end

function vertical = addVerticalLoad(vertical, name, magnitude, direction, x, z)
    arguments
        vertical
        name string
        magnitude double
        direction string {mustBeMember(direction,["down","up"])}
        x double
        z double
    end

    load = struct();
    load.name = char(name);
    load.value = magnitude;
    load.sign = strcmp(direction,"down") * 2 - 1; % +1 for down, -1 for up
    load.x = x;
    load.z = z;

    vertical(end+1) = load; %#ok<AGROW>
end

function zbar = centroidZTrapezoid(height, width_base, width_top)
    % Returns centroid measured from the base of the trapezoid
    if abs(width_base - width_top) < 1e-6
        zbar = height / 2;
        return;
    end

    a = width_base;
    b = width_top;
    zbar = height * (2*a + b) / (3 * (a + b));
end

function val = kPa(value)
    val = value; % 1 kN/m^2 equals 1 kPa
end
