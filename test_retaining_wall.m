%% Quick Test Script for Retaining Wall Design
% This script tests the complete design function with default parameters

clear; clc; close all;

fprintf('Testing retaining wall design scripts...\n\n');

% Test 1: Check if repmat works correctly
fprintf('Test 1: String repetition\n');
fprintf('%s\n', repmat('=', 1, 50));
fprintf('String repetition works correctly!\n');
fprintf('%s\n\n', repmat('=', 1, 50));

% Test 2: Run the complete design with minimal output
fprintf('Test 2: Running complete retaining wall design...\n\n');

try
    % Run the complete design script
    run('retaining_wall_complete_design.m');
    fprintf('\n? Complete design script executed successfully!\n');
catch ME
    fprintf('\n? Error in complete design script:\n');
    fprintf('   %s\n', ME.message);
    fprintf('   Line: %d\n', ME.stack(1).line);
end

fprintf('\nAll tests completed!\n');
