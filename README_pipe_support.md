# Pipe Support Design Calculation - MATLAB Code

## Overview
This MATLAB script performs comprehensive structural design calculations for a single-post pipe support system carrying a 12-inch water-filled stainless steel pipe.

## Design Configuration

### Pipe Specifications
- **Pipe Size**: 12-inch Schedule 10
- **Material**: Stainless Steel 304L
- **Outer Diameter**: 12.75 inches
- **Wall Thickness**: 0.165 inches
- **Content**: Water (full)

### Support Post
- **Section**: HSS 3" × 3" × 0.25" (square tube)
- **Height**: 1000 mm
- **Material**: ASTM A500 Grade B
- **Yield Strength**: 317 MPa

### Top Connection Plate
- **Dimensions**: 230 mm (L) × 100 mm (H) × 6 mm (thick)
- **Shape**: Curved to receive 12" pipe
- **Bolts**: 2× 7/8" diameter bolts (Grade A325)
- **Weld**: 6 mm fillet weld to post

### Base Plate & Anchors
- **Base Plate**: 8" × 8" × 0.25"
- **Anchor Bolts**: 4× Hilti 0.5" drop-in anchors
- **Concrete**: f'c = 25 MPa

## Calculations Performed

The script performs the following structural analysis:

1. **Section Properties**
   - Cross-sectional area, moment of inertia, section modulus, radius of gyration

2. **Load Calculations**
   - Pipe self-weight
   - Water weight
   - Support component weights
   - Factored loads (dead load factor = 1.25)
   - Lateral loads (assumed 5% of vertical)

3. **Post Structural Analysis**
   - Axial compression stress
   - Slenderness ratio check
   - Bending stress due to lateral load
   - Combined stress interaction
   - Lateral deflection check

4. **Weld Design**
   - Effective throat thickness
   - Axial and bending stress in weld
   - Combined stress check

5. **Top Bolt Design**
   - Shear stress from vertical load
   - Tension stress from overturning moment
   - Interaction check

6. **Top Plate Design**
   - Bearing stress under pipe load

7. **Base Plate & Anchor Design**
   - Concrete bearing stress
   - Anchor bolt shear stress
   - Anchor bolt tension (overturning)
   - Shear-tension interaction
   - Base plate bending stress

## Usage

### Running the Script
```matlab
% In MATLAB command window:
pipe_support_design
```

### Modifying Parameters
Edit the INPUT PARAMETERS section in the script to adjust:
- Pipe size and schedule
- Post dimensions
- Connection details
- Material properties
- Load factors

### Key Parameters You Can Adjust
```matlab
pipe_OD = 12.75;              % Pipe outside diameter (inches)
post_size = 3.0;              % Post size (inches)
post_wall_thick = 0.25;       % Post wall thickness (inches)
post_height = 1000;           % Height from floor (mm)
plate_thickness = 6;          % Top plate thickness (mm)
anchor_dia = 0.5;             % Anchor diameter (inches)
pipe_length_supported = 1000; % Span length (mm)
```

## Output

The script generates a detailed report including:

- Section properties for all components
- Load breakdown (pipe, water, self-weight)
- Stress utilization ratios for each component
- Deflection calculations
- Pass/fail status for each check
- Overall design status

### Sample Output
```
===============================================
   PIPE SUPPORT DESIGN CALCULATION
===============================================

1. SECTION PROPERTIES
   Post (HSS 3x3x0.25):
   - Cross-sectional area: 1419.35 mm²
   - Moment of inertia: 73.94 × 10⁴ mm⁴
   ...

2. LOAD CALCULATIONS
   Total pipe weight: 531.46 N (54.16 kg)
   Total water weight: 3157.87 N (321.99 kg)
   ...

DESIGN SUMMARY
   - Post axial stress: XX% utilized
   - Weld stress: XX% utilized
   ...

OVERALL DESIGN STATUS: ACCEPTABLE ✓
```

## Design Assumptions

1. **Loading**:
   - Dead load only (water-filled pipe)
   - Lateral load = 5% of vertical load (modify as needed)
   - No impact or dynamic loads

2. **Boundary Conditions**:
   - Post considered as fixed at base, free at top
   - Effective length factor K = 1.0

3. **Material Properties**:
   - Steel: E = 200,000 MPa
   - Concrete: f'c = 25 MPa

4. **Design Standards**:
   - Based on AISC/CSA steel design principles
   - ACI 318 for concrete bearing

## Safety Considerations

⚠️ **Important Notes**:
- This is a preliminary design calculation tool
- Verify all assumptions for your specific application
- Consider local building codes and standards
- Account for dynamic loads, impact factors, and seismic requirements
- Professional engineer review recommended before construction
- Field conditions may vary - site inspection required

## Modifications & Customization

### To add seismic or wind loads:
Modify the lateral load calculation section:
```matlab
H_lateral = 0.05 * P_dead;  % Change multiplier or add specific calculation
```

### To change design codes:
Update material strengths and allowable stress factors in the appropriate sections.

### To analyze different pipe schedules:
Update the pipe wall thickness parameter:
```matlab
pipe_wall_thick = 0.165;  % For SCH 10
% pipe_wall_thick = 0.250;  % For SCH 20
% pipe_wall_thick = 0.375;  % For SCH 40
```

## Requirements

- MATLAB R2016b or later
- No additional toolboxes required

## Contact & Support

For questions or modifications, refer to the inline comments in the code.
Each calculation section is clearly documented with formulas and assumptions.

---
**Version**: 1.0  
**Date**: November 6, 2025  
**Author**: Auto-generated design tool
