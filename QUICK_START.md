# Quick Start Guide - Retaining Wall Design

## ? All Scripts Fixed and Ready to Use!

The MATLAB syntax errors have been corrected. All scripts now use proper MATLAB syntax for string operations.

## ?? Start Here

### Option 1: Complete Design (Recommended)
```matlab
% Open and run this file:
retaining_wall_complete_design.m
```

### Option 2: Stability Only
```matlab
% Open and run this file:
retaining_wall_stability.m
```

### Option 3: Structural Design Only
```matlab
% Open and run this file:
retaining_wall_structural_design.m
```

## ?? Basic Usage Example

1. **Open** `retaining_wall_complete_design.m` in MATLAB
2. **Edit** the input parameters (starting at line 13):

```matlab
% Example: 6m high wall
H = 6.0;              % Wall height (m)
B = 4.0;              % Base width (m)
t_stem_base = 0.4;    % Stem thickness (m)
t_base = 0.6;         % Base thickness (m)
fc_prime = 30;        % Concrete: 30 MPa
fy = 400;             % Steel: 400 MPa
gamma_soil = 18.0;    % Soil density (kN/m?)
phi = 30;             % Soil friction angle (degrees)
```

3. **Run** the script: Press **F5** or click **Run**
4. **Review** the console output and figures

## ?? What You'll Get

### Console Output:
- ? Stability analysis (sliding, overturning, bearing)
- ? Structural design calculations
- ? Reinforcement requirements
- ? Pass/Fail indicators for all checks

### Figures (6 Drawings):
1. Wall section with dimensions
2. Stem reinforcement layout
3. Base slab reinforcement layout
4. Moment diagram
5. Shear diagram  
6. Bearing pressure distribution

## ?? Typical Design Workflow

### Step 1: Preliminary Sizing
Use these rules of thumb:
```
Base width:      B = 0.5 to 0.7 ? H
Stem thickness:  t = H / 16 to H / 12
Base thickness:  t = B / 8 to B / 6
```

### Step 2: Run Design
Execute the script with your preliminary dimensions.

### Step 3: Check Results
Look for:
- ? All safety factors ? required values
- ? Resultant within middle third
- ? Moment/shear capacity > demand

### Step 4: Optimize
If over-designed (< 70% utilization):
- Reduce section sizes
- Use smaller bar spacing

If under-designed:
- Increase dimensions
- Add more reinforcement

## ?? Common Adjustments

### Increase Overturning Safety:
```matlab
B = 4.5;          % Increase base width
t_heel = 2.8;     % Extend heel
```

### Increase Sliding Safety:
```matlab
B = 4.5;          % Increase base width
% Or add shear key in design
```

### Reduce Bearing Pressure:
```matlab
B = 5.0;          % Widen base
```

### Increase Moment Capacity:
```matlab
t_stem_base = 0.5;  % Thicken stem
% Or use higher strength materials
```

## ?? Design Parameters Reference

### Typical Soil Properties

| Soil Type | ? (kN/m?) | ? (deg) | Recommended Use |
|-----------|-----------|---------|-----------------|
| Loose sand | 16-18 | 28-30 | Not recommended |
| Dense sand | 18-20 | 32-36 | Good backfill |
| Sandy gravel | 18-20 | 34-40 | Excellent backfill |
| Silty sand | 17-19 | 28-32 | Fair backfill |

### Concrete Grades (Canadian)

| Grade | f'c (MPa) | Typical Use |
|-------|-----------|-------------|
| 25 MPa | 25 | Residential, light duty |
| 30 MPa | 30 | **Standard** (recommended) |
| 35 MPa | 35 | Heavy duty, high loads |

### Steel Grades (Canadian)

| Grade | fy (MPa) | Common Sizes |
|-------|----------|--------------|
| 400W | 400 | **Standard** (10M-35M) |
| 500W | 500 | High strength applications |

## ?? Important Safety Factors

The scripts check against these minimum values:

| Check | Minimum FS | CSA Standard |
|-------|------------|--------------|
| Overturning | 2.0 | Required |
| Sliding | 1.5 | Required |
| Bearing | 3.0 | Required |

## ?? Test Your Installation

Run the test script to verify everything works:

```matlab
test_retaining_wall
```

This will:
- Check MATLAB syntax compatibility
- Run a complete design example
- Verify all calculations
- Generate all figures

## ?? Additional Resources

- **Full documentation:** `README_retaining_wall_design.md`
- **Formula reference:** `DESIGN_FORMULAS.md`
- **CSA A23.3-19:** Design of Concrete Structures

## ?? Tips

1. **Start conservative:** Use larger dimensions initially, then optimize
2. **Check all outputs:** Don't just look at pass/fail, review actual values
3. **Use standard spacings:** 100, 150, 200, 250, 300 mm for easier construction
4. **Consider construction:** Practical bar arrangements matter
5. **Verify soil data:** Use actual geotechnical investigation results

## ?? Troubleshooting

### Error: "Ka is negative or complex"
**Cause:** Unrealistic soil parameters  
**Fix:** Check ? and ? values (? should be ? ?)

### Warning: "Shear reinforcement REQUIRED"
**Cause:** Section too thin  
**Fix:** Increase stem or base thickness by 100-150mm

### Error: "Cannot find suitable reinforcement"
**Cause:** Very high moment demand  
**Fix:** Increase section thickness or use higher strength materials

### Warning: "Resultant outside middle third"
**Cause:** Unstable bearing pressure distribution  
**Fix:** Increase base width or extend heel length

## ?? Support

For detailed help:
1. Check `README_retaining_wall_design.md` for comprehensive documentation
2. Review `DESIGN_FORMULAS.md` for calculation details
3. Consult CSA A23.3-19 for code requirements

---

## ? What Was Fixed

**Original Error:**
```matlab
fprintf('='*60);  % Python syntax - doesn't work in MATLAB
```

**Fixed Version:**
```matlab
fprintf('%s\n', repmat('=', 1, 60));  % Proper MATLAB syntax
```

All three main scripts have been corrected and are ready to use!

---

**Version:** 1.0 (Fixed)  
**Date:** October 2025  
**Status:** ? All scripts tested and working
