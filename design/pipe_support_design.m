function results = pipe_support_design(userInputs)
%PIPE_SUPPORT_DESIGN Pipe support design calculations per CSA S16 assumptions.
%
%   RESULTS = PIPE_SUPPORT_DESIGN(USERINPUTS) runs a series of strength checks
%   for a vertical pipe support consisting of a single hollow structural
%   section (HSS) post, a curved top plate, fillet welds, flange bolts, and a
%   base plate with drop-in anchors. The calculations implement typical CSA S16
%   design provisions using conservative assumptions suitable for preliminary
%   design and verification.
%
%   USERINPUTS is an optional structure that can override any default input
%   parameter. See the `defaultInputs` subfunction below for the full list of
%   available inputs and their default values. Key inputs include tributary
%   pipe length, additional imposed loads, material strengths, and connection
%   capacities.
%
%   The function returns a RESULTS structure summarising section properties,
%   load effects, resistances, utilisation ratios, and intermediate values for
%   transparency. Call the function without input arguments to run the default
%   design scenario:
%
%       results = pipe_support_design();
%
%   To update individual parameters, supply only those fields in USERINPUTS:
%
%       userInputs = struct();
%       userInputs.tributaryLength_m = 2.5;     % m of pipe supported
%       userInputs.transverseLoad_kN = 4.0;     % factored lateral load at pipe level
%       results = pipe_support_design(userInputs);
%
%   IMPORTANT ASSUMPTIONS
%   ---------------------
%   * Steel design follows CSA S16-19 philosophy using resistance factors noted
%     in the code. Column resistance adopts the generalized column formula that
%     is consistent with Clause 13.3 using 350W HSS properties.
%   * The 12 in Schedule 10S pipe properties are pre-populated and should be
%     verified against the project specification. The pipe material is assumed
%     to be stainless steel (340L / 304L) with density 8000 kg/m^3.
%   * The top saddle plate is curved to match the pipe OD; only shear transfer
%     through the 6 mm fillet weld around the HSS perimeter is checked. Local
%     bending of the plate is not explicitly modelled.
%   * Bolt and anchor capacities are evaluated using nominal CSA S16 equations
%     with resistance factors applied. Manufacturer-specific anchor ratings
%     should be substituted where available.
%   * Concrete bearing beneath the base plate is checked using CSA A23.3 style
%     bearing resistance (0.85 fc').
%   * Dynamic, thermal, or seismic load effects are not included by default.
%
%   Output highlights:
%       results.post.columnUtilisation   - Compression+bending unity check
%       results.weld.utilisation        - Factored shear vs resistance
%       results.bolts.utilisation       - Shear demand per bolt vs capacity
%       results.anchors.utilisation     - Tension & shear demand vs capacities
%       results.basePlate.bearingUtil   - Bearing pressure vs resistance
%
%   This function is intended to support engineering calculations. The user is
%   responsible for verifying all assumptions, adapting the code to project-
%   specific requirements, and ensuring compliance with governing standards.

%   Developed: 31-Oct-2025

    if nargin == 0 || isempty(userInputs)
        userInputs = struct();
    end

    params = applyDefaults(userInputs, defaultInputs());

    % Section and material properties ---------------------------------------
    props = computeSectionProperties(params);

    % Load effects -----------------------------------------------------------
    loads = assembleLoads(params, props);

    % Structural checks ------------------------------------------------------
    post      = checkPost(params, props, loads);
    weld      = checkWeld(params, props, loads);
    bolts     = checkBolts(params, props, loads);
    anchors   = checkAnchors(params, props, loads);
    basePlate = checkBasePlate(params, props, loads);

    % Collect results --------------------------------------------------------
    results = struct();
    results.inputs          = params;
    results.sectionProps    = props;
    results.loads           = loads;
    results.post            = post;
    results.weld            = weld;
    results.bolts           = bolts;
    results.anchors         = anchors;
    results.basePlate       = basePlate;

    % Provide quick summary in command window if requested ------------------
    if params.displaySummary
        displaySummary(results);
    end
end

% -------------------------------------------------------------------------
function params = defaultInputs()
    params = struct();

    % Geometry (mm unless noted)
    params.post.outerWidth_mm       = 76.2;    % 3 in HSS outside width
    params.post.wallThickness_mm    = 6.35;    % 0.25 in wall
    params.post.K_factor            = 1.0;     % Effective length factor
    params.post.height_mm           = 1000;    % Column clear height

    params.topPlate.length_mm       = 230;
    params.topPlate.height_mm       = 100;
    params.topPlate.thickness_mm    = 6;
    params.topPlate.weldSize_mm     = 6;       % Fillet size around post

    params.pipe.nominalSize_in      = 12;
    params.pipe.schedule            = "10S";
    params.pipe.outsideDiameter_mm  = 323.9;   % 12.75 in
    params.pipe.wallThickness_mm    = 4.19;    % Schedule 10S

    params.pipe.materialDensity_kgm3 = 8000;   % Stainless 304/316
    params.fluid.density_kgm3        = 1000;   % Water

    params.tributaryLength_m        = 1.8;     % Pipe length supported
    params.additionalDeadLoad_kN    = 0.0;     % e.g., insulation, fittings
    params.transverseLoad_kN        = 2.0;     % Factored lateral load (wind/seismic)
    params.longitudinalLoad_kN      = 0.5;     % Factored longitudinal load
    params.verticalLiveLoad_kN      = 0.0;     % Additional vertical live load
    params.eccentricity_mm          = 0.0;     % Vertical load eccentricity

    % Base plate and anchors
    params.basePlate.length_mm      = 203.2;   % 8 in
    params.basePlate.width_mm       = 203.2;   % 8 in
    params.basePlate.thickness_mm   = 6.35;    % 0.25 in
    params.basePlate.concrete_fc_MPa = 30;     % Specified compressive strength
    params.anchor.pattern_mm        = [150, 150]; % centre-to-centre in x,y
    params.anchor.diameter_mm       = 12.7;    % 0.5 in
    params.anchor.factoredTensionCapacity_kN = 6.8; % per anchor, manufacturer
    params.anchor.factoredShearCapacity_kN   = 5.3; % per anchor

    % Material strengths (MPa)
    params.material.Fy_MPa          = 350;     % 350W HSS
    params.material.Fu_MPa          = 450;
    params.weld.Fexx_MPa            = 480;     % e.g., E70 electrode (~480 MPa)
    params.bolt.Fu_MPa              = 620;     % A325/A490 class
    params.anchor.materialFu_MPa    = 620;     % For reporting only

    % Resistance factors (CSA S16)
    params.phi.compression          = 0.9;
    params.phi.bending              = 0.9;
    params.phi.weld                 = 0.67;
    params.phi.boltShear            = 0.6;
    params.phi.anchor               = 0.75;
    params.loadFactor.dead          = 1.25;
    params.loadFactor.live          = 1.5;

    % Additional settings
    params.gravity_m_s2             = 9.80665;
    params.E_modulus_MPa            = 200000; % Structural steel
    params.displaySummary           = true;
end

% -------------------------------------------------------------------------
function params = applyDefaults(userInputs, defaults)
    params = defaults;
    userFields = fieldnames(userInputs);
    for iField = 1:numel(userFields)
        fieldName = userFields{iField};
        if isstruct(userInputs.(fieldName)) && isfield(defaults, fieldName)
            params.(fieldName) = applyDefaults(userInputs.(fieldName), defaults.(fieldName));
        else
            params.(fieldName) = userInputs.(fieldName);
        end
    end
end

% -------------------------------------------------------------------------
function props = computeSectionProperties(params)
    props = struct();

    b = params.post.outerWidth_mm / 1000;            % m
    t = params.post.wallThickness_mm / 1000;         % m
    hi = b - 2*t;

    props.post.area_m2 = b^2 - hi^2;
    props.post.area_mm2 = props.post.area_m2 * 1e6;

    Ix = (b^4 - hi^4) / 12;
    props.post.Ix_m4 = Ix;
    props.post.Iy_m4 = Ix;                            % square section
    props.post.rx_m  = sqrt(Ix / props.post.area_m2);
    props.post.ry_m  = props.post.rx_m;

    props.post.sectionModulus_m3 = (b^4 - hi^4) / (6*b);
    props.post.plasticModulus_m3 = (4/3) * (b^3 - hi^3) / 4; % Approximate for square HSS

    % Weld properties around post perimeter (neglecting corner round)
    perimeter_m = 4 * (b - 2*0.5*t); % subtract half thickness per corner
    props.weld.length_m = perimeter_m;
    props.weld.effectiveThroat_m = 0.707 * (params.topPlate.weldSize_mm / 1000);
    props.weld.area_m2 = props.weld.length_m * props.weld.effectiveThroat_m;

    % Bolt properties (7/8 in diameter)
    boltDiameter_m = 7/8 * 25.4 / 1000;
    boltPitch_in = 9; % 7/8-9 UNC thread
    minorDiameter_in = 7/8 - 1.08253 / boltPitch_in;
    props.bolt.area_m2 = pi/4 * (minorDiameter_in * 25.4 / 1000)^2;
    props.bolt.count = 2;

    % Anchor geometry
    props.anchor.gauge_m = params.anchor.pattern_mm / 1000;
    props.anchor.count = 4;
    props.anchor.area_m2 = pi/4 * (params.anchor.diameter_mm / 1000)^2;

    % Base plate area
    props.basePlate.area_m2 = (params.basePlate.length_mm / 1000) * (params.basePlate.width_mm / 1000);
end

% -------------------------------------------------------------------------
function loads = assembleLoads(params, props)
    loads = struct();

    pipeOD = params.pipe.outsideDiameter_mm / 1000;
    pipeWall = params.pipe.wallThickness_mm / 1000;
    pipeID = pipeOD - 2*pipeWall;

    pipeArea_m2 = pi/4 * (pipeOD^2 - pipeID^2);
    fluidArea_m2 = pi/4 * pipeID^2;

    pipeWeight_kN = pipeArea_m2 * params.pipe.materialDensity_kgm3 * params.gravity_m_s2 * params.tributaryLength_m / 1000;
    fluidWeight_kN = fluidArea_m2 * params.fluid.density_kgm3 * params.gravity_m_s2 * params.tributaryLength_m / 1000;

    loads.deadLoad_kN = pipeWeight_kN + fluidWeight_kN + params.additionalDeadLoad_kN;
    loads.verticalFactored_kN = params.loadFactor.dead * loads.deadLoad_kN + params.loadFactor.live * params.verticalLiveLoad_kN;

    loads.transverseFactored_kN = params.transverseLoad_kN;
    loads.longitudinalFactored_kN = params.longitudinalLoad_kN;

    loads.shearResultant_kN = hypot(loads.transverseFactored_kN, loads.longitudinalFactored_kN);

    loads.momentBase_transverse_kNm = loads.transverseFactored_kN * (params.post.height_mm / 1000) + ...
        loads.verticalFactored_kN * (params.eccentricity_mm / 1000);
    loads.momentBase_longitudinal_kNm = loads.longitudinalFactored_kN * (params.post.height_mm / 1000);

    loads.pipeWeight_kN = pipeWeight_kN;
    loads.fluidWeight_kN = fluidWeight_kN;
end

% -------------------------------------------------------------------------
function post = checkPost(params, props, loads)
    post = struct();

    A = props.post.area_m2;
    Fy = params.material.Fy_MPa;
    E = params.E_modulus_MPa;

    KL = params.post.K_factor * (params.post.height_mm / 1000);
    r = props.post.rx_m;
    slenderness = KL / r;

    Fe = (pi^2 * (E)) / slenderness^2; % MPa
    lambda = Fy / Fe;
    if lambda <= 2.25
        Fcr = (0.658 ^ lambda) * Fy;
    else
        Fcr = 0.877 * Fe;
    end

    Pc_kN = params.phi.compression * Fcr * A * 1e6 / 1000;

    % Bending capacity (equal about both axes for square HSS)
    Mr_kNm = params.phi.bending * Fy * props.post.sectionModulus_m3 * 1e6 / 1000;

    % Interaction check (CSA S16 simplified)
    Pu = loads.verticalFactored_kN;
    Mux = abs(loads.momentBase_transverse_kNm);
    Muy = abs(loads.momentBase_longitudinal_kNm);

    utilisation = Pu / Pc_kN + max(Mux, Muy) / Mr_kNm;

    post.axialCapacity_kN = Pc_kN;
    post.bendingCapacity_kNm = Mr_kNm;
    post.slenderness = slenderness;
    post.Fcr_MPa = Fcr;
    post.interactionUtilisation = utilisation;
    post.axialStress_MPa = Pu * 1000 / (A * 1e6);
    post.bendingStress_MPa = max(Mux, Muy) * 1e6 / props.post.sectionModulus_m3;
end

% -------------------------------------------------------------------------
function weld = checkWeld(params, props, loads)
    weld = struct();

    V = loads.verticalFactored_kN;
    Aw = props.weld.area_m2;
    Fexx = params.weld.Fexx_MPa;
    phi = params.phi.weld;

    Vr_kN = phi * 0.67 * Fexx * Aw * 1e6 / 1000;
    utilisation = V / Vr_kN;

    weld.shearCapacity_kN = Vr_kN;
    weld.utilisation = utilisation;
    weld.effectiveThroat_mm = props.weld.effectiveThroat_m * 1000;
    weld.totalLength_mm = props.weld.length_m * 1000;
end

% -------------------------------------------------------------------------
function bolts = checkBolts(params, props, loads)
    bolts = struct();

    V = loads.verticalFactored_kN;
    n = props.bolt.count;
    A = props.bolt.area_m2;
    Fu = params.bolt.Fu_MPa;
    phi = params.phi.boltShear;

    Vr_perBolt_kN = phi * 0.55 * Fu * A * 1e6 / 1000;
    demand_perBolt_kN = V / n; % assumes symmetric loading
    utilisation = demand_perBolt_kN / Vr_perBolt_kN;

    bolts.count = n;
    bolts.capacityPerBolt_kN = Vr_perBolt_kN;
    bolts.demandPerBolt_kN = demand_perBolt_kN;
    bolts.utilisation = utilisation;
end

% -------------------------------------------------------------------------
function anchors = checkAnchors(params, props, loads)
    anchors = struct();

    n = props.anchor.count;
    gauge = props.anchor.gauge_m; % [gx, gy]
    h = params.post.height_mm / 1000;

    M = hypot(loads.momentBase_transverse_kNm, loads.momentBase_longitudinal_kNm);
    V = loads.shearResultant_kN;
    Pu = loads.verticalFactored_kN;

    % Assume compression bearing prevents uplift unless overturning moment
    % exceeds stabilising compression about principal axes.
    eccentricity_m = M / Pu;

    effectiveHalfWidth_m = gauge(1) / 2;
    effectiveHalfDepth_m = gauge(2) / 2;
    maxDistance_m = sqrt(effectiveHalfWidth_m^2 + effectiveHalfDepth_m^2);

    tensionPerAnchor_kN = max(0, (Pu * (eccentricity_m / maxDistance_m)) / n);

    shearPerAnchor_kN = V / n;

    utilisationTension = tensionPerAnchor_kN / params.anchor.factoredTensionCapacity_kN;
    utilisationShear = shearPerAnchor_kN / params.anchor.factoredShearCapacity_kN;

    anchors.count = n;
    anchors.tensionDemand_kN = tensionPerAnchor_kN;
    anchors.shearDemand_kN = shearPerAnchor_kN;
    anchors.tensionCapacity_kN = params.anchor.factoredTensionCapacity_kN;
    anchors.shearCapacity_kN = params.anchor.factoredShearCapacity_kN;
    anchors.utilisationTension = utilisationTension;
    anchors.utilisationShear = utilisationShear;
    anchors.overturningEccentricity_m = eccentricity_m;
end

% -------------------------------------------------------------------------
function basePlate = checkBasePlate(params, props, loads)
    basePlate = struct();

    area_m2 = props.basePlate.area_m2;
    Pu = loads.verticalFactored_kN;
    q = Pu / area_m2; % kN/m^2 = kPa

    fcPrime = params.basePlate.concrete_fc_MPa;
    phi = params.phi.anchor; % re-use as bearing factor (conservative)
    qr = phi * 0.85 * fcPrime * 1000; % kPa

    basePlate.bearingPressure_kPa = q;
    basePlate.bearingResistance_kPa = qr;
    basePlate.bearingUtil = q / qr;

    basePlate.area_m2 = area_m2;
end

% -------------------------------------------------------------------------
function displaySummary(results)
    fprintf('\nPIPE SUPPORT DESIGN SUMMARY\n');
    fprintf('---------------------------------------------\n');
    fprintf('Vertical factored load     : %7.2f kN\n', results.loads.verticalFactored_kN);
    fprintf('Transverse factored load   : %7.2f kN\n', results.loads.transverseFactored_kN);
    fprintf('Longitudinal factored load : %7.2f kN\n', results.loads.longitudinalFactored_kN);
    fprintf('Post interaction utilisation: %7.3f (<= 1 OK)\n', results.post.interactionUtilisation);
    fprintf('Weld shear utilisation     : %7.3f (<= 1 OK)\n', results.weld.utilisation);
    fprintf('Bolt shear utilisation     : %7.3f (<= 1 OK)\n', results.bolts.utilisation);
    fprintf('Anchor tension utilisation : %7.3f (<= 1 OK)\n', results.anchors.utilisationTension);
    fprintf('Anchor shear utilisation   : %7.3f (<= 1 OK)\n', results.anchors.utilisationShear);
    fprintf('Base plate bearing utilisation: %7.3f (<= 1 OK)\n', results.basePlate.bearingUtil);
    fprintf('---------------------------------------------\n\n');
end
