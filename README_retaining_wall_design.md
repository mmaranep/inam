# Retaining Wall Design Suite - Canadian Codes (CSA A23.3-19)

## Overview

This suite provides comprehensive MATLAB tools for designing cantilever retaining walls according to Canadian standards (CSA A23.3-19 - Design of Concrete Structures).

## Files Included

### 1. `retaining_wall_stability.m`
**Purpose:** Global stability analysis only
- Sliding stability check
- Overturning stability check
- Bearing capacity analysis
- Middle-third rule verification
- Detailed pressure distribution visualization

**Use this when:** You need quick stability verification or preliminary design.

### 2. `retaining_wall_structural_design.m`
**Purpose:** Complete structural design per CSA A23.3-19
- Stem flexural and shear design
- Toe slab design
- Heel slab design
- Development length checks
- Crack control verification
- Temperature and shrinkage reinforcement
- Detailed reinforcement drawings

**Use this when:** You need full structural design with reinforcement details.

### 3. `retaining_wall_complete_design.m` ? **RECOMMENDED**
**Purpose:** Combined stability + structural design
- All features from both above scripts
- Comprehensive design report
- Professional-quality drawings
- Complete reinforcement summary
- Material specifications

**Use this when:** You need a complete design package (most common use case).

## Quick Start

### Basic Usage

1. Open `retaining_wall_complete_design.m` in MATLAB
2. Modify input parameters (lines 13-55):
```matlab
% Geometry
H = 6.0;              % Wall height (m)
B = 4.0;              % Base width (m)
t_stem_base = 0.4;    % Stem thickness (m)
t_base = 0.6;         % Base thickness (m)

% Soil properties
gamma_soil = 18.0;    % Backfill unit weight (kN/m?)
phi = 30;             % Friction angle (degrees)

% Concrete properties
fc_prime = 30;        % Concrete strength (MPa)
fy = 400;             % Steel yield strength (MPa)
```

3. Run the script: Press **F5** or click **Run**
4. Review console output and generated figures

## Input Parameters Guide

### Geometry Parameters

| Parameter | Description | Typical Range | Units |
|-----------|-------------|---------------|-------|
| `H` | Total wall height | 3.0 - 10.0 | m |
| `B` | Base width | H/2 to 2H/3 | m |
| `t_stem_base` | Stem thickness at base | 0.3 - 0.6 | m |
| `t_base` | Base slab thickness | 0.4 - 0.8 | m |
| `t_toe` | Toe length | B/3 to B/2 | m |
| `cover` | Concrete cover | 0.075 | m |

### Soil Parameters

| Parameter | Description | Typical Values | Units |
|-----------|-------------|----------------|-------|
| `gamma_soil` | Backfill unit weight | 16-20 (sand/gravel) | kN/m? |
| `phi` | Internal friction angle | 28-38? (granular) | degrees |
| `c` | Cohesion | 0 (granular), 10-50 (cohesive) | kPa |
| `q` | Surcharge load | 0-20 | kPa |

### Material Parameters

| Parameter | Description | Typical Values | Units |
|-----------|-------------|----------------|-------|
| `fc_prime` | Concrete strength | 25, 30, 35, 40 | MPa |
| `fy` | Steel yield strength | 400, 500 | MPa |
| `gamma_concrete` | Concrete unit weight | 24 | kN/m? |

### Safety Factors (CSA Standards)

| Check | Minimum FS |
|-------|------------|
| Sliding | 1.5 |
| Overturning | 2.0 |
| Bearing Capacity | 3.0 |

## Output Interpretation

### Console Output Structure

#### Part A: Global Stability
```
OVERTURNING CHECK:
  FS = 2.45 (Required: 2.0) [PASS] ?

SLIDING CHECK:
  FS = 1.68 (Required: 1.5) [PASS] ?

BEARING CAPACITY CHECK:
  FS = 3.52 (Required: 3.0) [PASS] ?
```

#### Part B: Structural Design
```
STEM DESIGN:
  PROVIDED: 20M @ 200 mm c/c
  Mr = 125.5 kN?m/m (87.2% utilized)

TOE SLAB:
  PROVIDED: 15M @ 200 mm c/c (bottom)
  
HEEL SLAB:
  PROVIDED: 20M @ 150 mm c/c (top)
```

### Visual Outputs

The script generates 6 key drawings:

1. **Wall Section** - Overall geometry with dimensions
2. **Stem Reinforcement** - Bar layout in stem
3. **Base Reinforcement** - Toe and heel bar arrangement
4. **Moment Diagram** - Bending moment distribution
5. **Shear Diagram** - Shear force distribution
6. **Bearing Pressure** - Soil pressure distribution

## Design Workflow

### Step 1: Preliminary Sizing

**Rule of thumb for initial dimensions:**
- Base width: `B = 0.5H to 0.7H`
- Stem base thickness: `t_stem = H/24 to H/16`
- Base thickness: `t_base = B/8 to B/6`
- Toe length: `t_toe = B/3 to B/2`

### Step 2: Stability Check

Run the script and check:
- ? All stability factors of safety met
- ? Resultant within middle third
- ? Bearing pressure < allowable

If any check fails, adjust:
- **Overturning fails:** Increase `B` or `t_heel`
- **Sliding fails:** Increase `B`, add key, or use passive resistance
- **Bearing fails:** Increase `B` or improve foundation soil

### Step 3: Structural Design

Check:
- ? Moment capacity > factored moment
- ? Shear capacity > factored shear
- ? Development length adequate
- ? Crack control satisfied

If any check fails:
- Increase member thickness
- Add more reinforcement
- Use higher grade steel/concrete

### Step 4: Optimization

Fine-tune the design for economy:
- Reduce over-designed sections (< 70% utilization)
- Standardize reinforcement spacing (100, 150, 200, 250 mm)
- Minimize number of different bar sizes

## Design Checks Summary

### CSA A23.3-19 Code Compliance

The scripts implement:

? **Clause 10.5** - Minimum flexural reinforcement
? **Clause 11.3** - Shear resistance of members without shear reinforcement  
? **Clause 12.2** - Development of reinforcement
? **Clause 10.6** - Crack control
? **Clause 7.8** - Temperature and shrinkage reinforcement
? **Material resistance factors:** ?c = 0.65, ?s = 0.85
? **Load factors:** ?D = 1.25, ?L = 1.5

## Common Issues & Solutions

### Issue 1: "Shear reinforcement REQUIRED"
**Solution:** Increase stem or base thickness by 100-150mm

### Issue 2: "Development length INADEQUATE"
**Solutions:**
- Use 90? or 180? hooks at bar ends
- Increase cover/embedment length
- Use smaller diameter bars

### Issue 3: "Resultant outside middle third"
**Solutions:**
- Increase base width `B`
- Extend heel length
- Reduce toe length
- Add backfill weight on heel

### Issue 4: Sliding failure
**Solutions:**
- Increase base width
- Add shear key under base slab
- Improve foundation soil
- Consider passive resistance (if reliable)

## Design Example

### Example: 6m High Retaining Wall

**Given:**
- Wall height: 6.0 m
- Backfill: Sandy gravel (? = 18 kN/m?, ? = 30?)
- Surcharge: 10 kPa
- Foundation soil: ?f = 28?
- Concrete: f'c = 30 MPa
- Steel: fy = 400 MPa

**Results:**
- Base width: 4.0 m
- Stem thickness: 400 mm
- Base thickness: 600 mm
- **Reinforcement:**
  - Stem: 20M @ 200 mm c/c (vertical)
  - Toe: 15M @ 200 mm c/c (bottom)
  - Heel: 20M @ 150 mm c/c (top)
  - Temperature: 10M @ 300 mm c/c

**Safety Factors:**
- Overturning: FS = 2.45 ?
- Sliding: FS = 1.68 ?
- Bearing: FS = 3.52 ?

## Advanced Features

### Custom Bar Selection

Modify `bar_sizes` and `bar_areas` arrays to match available bar stock:
```matlab
bar_sizes = [10, 15, 20, 25, 30, 35];  % mm
bar_areas = [100, 200, 300, 500, 700, 1000];  % mm?
```

### Load Factor Modifications

For different load combinations (CSA A23.3 Cl. 8.3):
```matlab
% Principal load combination
alpha_D = 1.25;
alpha_L = 1.5;

% Companion load combination
alpha_D = 1.25;
alpha_L = 1.0;
```

### Seismic Considerations

For seismic regions, increase factors:
```matlab
FS_sliding_min = 1.8;      % Instead of 1.5
FS_overturning_min = 2.5;  % Instead of 2.0
```

## Bar Designation (Canadian Practice)

| Bar Size | Diameter (mm) | Area (mm?) | Common Use |
|----------|---------------|------------|------------|
| 10M | 11.3 | 100 | Temperature/shrinkage |
| 15M | 16.0 | 200 | Light reinforcement |
| 20M | 19.5 | 300 | Main reinforcement |
| 25M | 25.2 | 500 | Heavy reinforcement |
| 30M | 29.9 | 700 | Heavy loads |
| 35M | 35.7 | 1000 | Very heavy loads |

## References

1. **CSA A23.3-19:** Design of Concrete Structures
2. **CSA A23.1-19:** Concrete Materials and Methods of Concrete Construction
3. **NBCC 2020:** National Building Code of Canada
4. **Foundation Engineering Handbook** - Fang & Winterkorn

## Support & Modifications

### To modify for different codes:

**US (ACI 318):**
- Change ?c = 0.65 ? 0.75
- Change ?s = 0.85 ? 0.90
- Update load factors per ACI 318

**Eurocode (EN 1992):**
- Implement partial safety factors
- Use different material properties
- Update reinforcement designations

## License & Disclaimer

These scripts are for educational and professional use. Always verify designs with:
- Licensed professional engineer
- Local building codes
- Geotechnical investigation report
- Peer review for critical structures

---

**Version:** 1.0  
**Last Updated:** October 2025  
**Compatible with:** MATLAB R2018a and later

For questions or suggestions, please contact your structural engineering team.
