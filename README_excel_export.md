# Pipe Support Design - Excel Export Tool

## Overview
This MATLAB script exports the complete pipe support design calculation to a comprehensive Excel spreadsheet with multiple formatted sheets.

## Output File
**Filename**: `Pipe_Support_Design_Calculation.xlsx`

The Excel file contains **9 detailed sheets**:

### Sheet 1: Input Parameters
- Complete list of all design inputs
- Pipe properties (12" SCH 10 SS304L)
- Post specifications (HSS 3×3×0.25)
- Top plate dimensions
- Bolt specifications
- Weld details
- Base plate and anchor properties
- Material properties
- Load factors

### Sheet 2: Section Properties
- Cross-sectional area
- Moment of inertia
- Section modulus
- Radius of gyration
- All geometric properties of the post

### Sheet 3: Load Calculations
- Pipe weight (empty)
- Water weight (full)
- Component self-weights
- Total dead load
- Factored loads
- Lateral load assumptions
- Load summary

### Sheet 4: Post Analysis
- Axial compression stress check
- Slenderness ratio verification
- Bending stress analysis
- Combined stress interaction
- Deflection calculations
- Status indicators (OK/FAIL) for each check

### Sheet 5: Weld Check
- Weld properties (6mm fillet)
- Throat thickness
- Effective weld area
- Axial stress in weld
- Moment stress in weld
- Combined stress check
- Utilization ratios

### Sheet 6: Top Bolt Check
- Bolt properties (7/8" diameter)
- Shear stress analysis
- Tension stress from moment
- Allowable stresses
- Utilization percentages
- Pass/fail indicators

### Sheet 7: Top Plate Check
- Plate dimensions
- Bearing stress under pipe load
- Allowable bearing stress
- Utilization ratio
- Status check

### Sheet 8: Base Plate & Anchors
- Concrete bearing stress
- Anchor bolt shear analysis
- Anchor bolt tension (overturning)
- Shear-tension interaction
- Base plate bending stress
- Complete anchor design verification

### Sheet 9: Summary
- Project information
- Loading summary
- Complete utilization table for all components
- Overall design status (ACCEPTABLE/REQUIRES REVISION)
- Design notes and recommendations
- Professional review reminders

## Usage

### Running the Script

1. **Open MATLAB** and navigate to the workspace directory:
```matlab
cd /path/to/workspace
```

2. **Run the export script**:
```matlab
pipe_support_to_excel
```

3. **Check the output**:
The script will create `Pipe_Support_Design_Calculation.xlsx` in the current directory.

### Sample Output
```
===============================================
   PIPE SUPPORT DESIGN - EXCEL EXPORT
===============================================

Creating Excel file: Pipe_Support_Design_Calculation.xlsx

Writing Sheet 1: Input Parameters...
Writing Sheet 2: Section Properties...
Writing Sheet 3: Load Calculations...
Writing Sheet 4: Post Analysis...
Writing Sheet 5: Weld Check...
Writing Sheet 6: Top Bolt Check...
Writing Sheet 7: Top Plate Check...
Writing Sheet 8: Base Plate & Anchors...
Writing Sheet 9: Summary...

===============================================
   Excel file created successfully!
   Filename: Pipe_Support_Design_Calculation.xlsx
===============================================

Overall Design Status: ACCEPTABLE ✓
```

## Features

### Automatic Status Indicators
Each check includes automatic status indicators:
- ✓ **OK** - Component passes the design check
- ✗ **FAIL** - Component requires revision

### Utilization Ratios
All components show utilization percentages:
- < 100% = Acceptable
- ≥ 100% = Over-stressed (requires revision)

### Comprehensive Documentation
Each sheet includes:
- Parameter symbols
- Units
- Descriptions
- Formulas (where applicable)
- Status checks

## Modifying Parameters

To modify design parameters, edit the **INPUT PARAMETERS** section at the top of the script:

```matlab
%% INPUT PARAMETERS

% Pipe Properties
pipe_OD = 12.75;              % Change pipe size
pipe_schedule = 10;           % Change schedule
pipe_length_supported = 1000; % Change span length

% Post Properties
post_size = 3.0;              % Change post size
post_height = 1000;           % Change height

% Plate Properties
plate_thickness = 6;          % Change plate thickness

% Anchor Properties
anchor_dia = 0.5;             % Change anchor size
anchor_num = 4;               % Change number of anchors

% Load Factors
SF_dead = 1.25;              % Adjust safety factor
```

After modifying parameters, simply re-run the script to generate a new Excel file with updated calculations.

## Excel File Structure

### Data Organization
- **Column A**: Parameter names/descriptions
- **Column B**: Symbols (engineering notation)
- **Column C**: Values
- **Column D**: Units
- **Column E**: Status/Notes

### Status Indicators
The script automatically adds status indicators:
- "OK" = Design check passed
- "FAIL" = Design check failed
- "PASS" = Alternative pass indicator

## Key Calculations Exported

### Post Design
- Axial stress vs. 0.85×Fy
- Bending stress vs. 0.66×Fy
- Combined stress interaction ≤ 1.0
- Slenderness KL/r ≤ 200
- Deflection ≤ H/200

### Connection Design
- Weld stress vs. 0.3×Fe
- Bolt shear vs. 0.4×Fu
- Bolt tension vs. 0.75×Fu
- Plate bearing vs. 0.9×Fy

### Foundation Design
- Concrete bearing vs. 0.35×f'c
- Anchor shear and tension
- Interaction: (V/Va)² + (T/Ta)² ≤ 1.0
- Base plate bending vs. 0.75×Fy

## Professional Use

### Documentation Package
The Excel file serves as complete documentation for:
- Design submittals
- Building permit applications
- Engineering review
- Construction reference
- Project records

### Review Checklist
Before finalizing the design, verify:
- [ ] All status indicators show "OK"
- [ ] Overall design status is "ACCEPTABLE"
- [ ] Load assumptions are appropriate
- [ ] Material grades are correct
- [ ] Anchor embedment is adequate
- [ ] Edge distances meet requirements
- [ ] Professional engineer has reviewed

## Advantages of Excel Format

1. **Portability**: Share with contractors, reviewers, inspectors
2. **Editing**: Can add notes, comments, or modifications
3. **Printing**: Easy to print for field use
4. **Integration**: Import into other software or reports
5. **Archiving**: Standard format for long-term storage
6. **Verification**: Easy for others to check calculations

## Troubleshooting

### File Already Open Error
If you get an error about the file being in use:
1. Close the Excel file
2. Re-run the script

The script automatically deletes the old file before creating a new one.

### Missing Excel
If MATLAB cannot write to Excel:
- Ensure Microsoft Excel is installed (Windows)
- Or use LibreOffice Calc (Linux/Mac)
- MATLAB's `writecell` function works with both

### Large File Size
If the Excel file is very large:
- This is normal for detailed calculations
- Typical size: 100-500 KB
- Compress before emailing if needed

## Additional Scripts

### To run both calculation and export:
```matlab
% Run calculation with console output
pipe_support_design

% Export to Excel
pipe_support_to_excel
```

### Batch Processing
To analyze multiple configurations:
```matlab
% Loop through different parameters
for height = [800, 1000, 1200]
    post_height = height;
    pipe_support_to_excel;
    movefile('Pipe_Support_Design_Calculation.xlsx', ...
             sprintf('Design_Height_%dmm.xlsx', height));
end
```

## Requirements

- MATLAB R2019a or later (for `writecell` function)
- Microsoft Excel, LibreOffice, or compatible spreadsheet software
- No additional toolboxes required

## Version History

- **v1.0** (Nov 6, 2025): Initial release
  - 9 comprehensive sheets
  - All structural checks included
  - Automatic status indicators

## Support

For questions or issues:
1. Check inline code comments
2. Verify input parameters
3. Review MATLAB console output
4. Check Excel file for error messages

---
**Note**: This Excel export is for documentation and verification purposes. The underlying calculations follow AISC/CSA standards, but professional engineering judgment and code compliance verification are required before construction.
