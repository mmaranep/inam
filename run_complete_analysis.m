%% COMPLETE PIPE SUPPORT ANALYSIS
% This script runs both the console calculation and Excel export
% Author: Auto-generated design tool
% Date: November 6, 2025

clear all;
clc;

fprintf('====================================================\n');
fprintf('   PIPE SUPPORT DESIGN - COMPLETE ANALYSIS\n');
fprintf('====================================================\n\n');

fprintf('This script will:\n');
fprintf('1. Run the design calculation with console output\n');
fprintf('2. Export all results to Excel spreadsheet\n\n');

% Prompt user
fprintf('Press any key to start...\n');
pause;

%% STEP 1: Run Console Calculation
fprintf('\n====================================================\n');
fprintf('STEP 1: Running Design Calculation...\n');
fprintf('====================================================\n\n');

try
    pipe_support_design;
    fprintf('\n✓ Design calculation completed successfully!\n\n');
catch ME
    fprintf('\n✗ Error in design calculation: %s\n\n', ME.message);
    return;
end

% Pause to let user review
fprintf('Press any key to continue to Excel export...\n');
pause;

%% STEP 2: Export to Excel
fprintf('\n====================================================\n');
fprintf('STEP 2: Exporting to Excel...\n');
fprintf('====================================================\n\n');

try
    pipe_support_to_excel;
    fprintf('\n✓ Excel export completed successfully!\n\n');
catch ME
    fprintf('\n✗ Error in Excel export: %s\n\n', ME.message);
    return;
end

%% COMPLETION MESSAGE
fprintf('====================================================\n');
fprintf('   ANALYSIS COMPLETE!\n');
fprintf('====================================================\n\n');

fprintf('Files created:\n');
fprintf('  • Pipe_Support_Design_Calculation.xlsx\n\n');

fprintf('Next steps:\n');
fprintf('  1. Open the Excel file to review detailed results\n');
fprintf('  2. Check the Summary sheet for overall status\n');
fprintf('  3. Verify all components show "OK" status\n');
fprintf('  4. Review assumptions and load cases\n');
fprintf('  5. Have a professional engineer review the design\n\n');

fprintf('Thank you for using the Pipe Support Design Tool!\n');
fprintf('====================================================\n');
