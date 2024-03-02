
function run_template_and_precision_convergence_tests()

flag_run_linear = 1;
flag_run_nearest = 1;
flag_run_precision = 1;

%%
% The known recording location of the electrode specified in the JSON file
% is 15 mm (15000 um). Meanwhile, the fiber spans from -70 mm to +30 mm. Thus,
% specify starting AP locations given a maximum of 85 mm conduction
% distance
KNOWN_CONDUCTION_DISTANCE_MM = 11; % [mm]

desired_conduction_distances_mm = [6,11,21,41,81]; % [mm]
% Define x_bounds for plotting for each conduction distance; if [], set
% automatically
x_bounds_myelinated = {[0,2],[],[0 10],[],[1.8,22.8]};
x_bounds_unmyelinated = {[0,22],[],[17,70],[],[75,190]};

% Run the fiber diameter precision convergence test across multiple different
% location_where_first_AP_occurred (i.e., different conduction distances)
% and plot the
parameter_sweep_values = KNOWN_CONDUCTION_DISTANCE_MM - desired_conduction_distances_mm;
parameter_sweep_variable_name = 'stim_location_change_mm';
title_labels = arrayfun(@(x) num2str(x,'%d mm'), round(desired_conduction_distances_mm),...
    'UniformOutput',false);

if (flag_run_nearest)
    % Run number of templates convergence tests
    % Nearest neighbor
    % myelinated
    parameter_sweep_xlabel = 'number of templates';
    reference_CNAP_label = '193';
    run_parameter_sweep_test('../JSON_input_params/RUN_CNAP_convergence_test_number_of_templates_nearest1.json',...
        parameter_sweep_values,parameter_sweep_variable_name,parameter_sweep_xlabel,title_labels,reference_CNAP_label,x_bounds_myelinated);

    reference_CNAP_label = '97';
    run_parameter_sweep_test('../JSON_input_params/RUN_CNAP_convergence_test_number_of_templates_nearest2.json',...
        parameter_sweep_values,parameter_sweep_variable_name,parameter_sweep_xlabel,title_labels,reference_CNAP_label,x_bounds_unmyelinated);
end

if (flag_run_precision)
    % Run precision convergence tests
    parameter_sweep_xlabel = 'precision (\mum)';
    reference_CNAP_label = '0.000001';
    %myel
    run_parameter_sweep_test('../JSON_input_params/RUN_CNAP_convergence_test_fiber_diameter_precision1.json',...
        parameter_sweep_values,parameter_sweep_variable_name,parameter_sweep_xlabel,title_labels,reference_CNAP_label,x_bounds_myelinated);
    %unmyel
    run_parameter_sweep_test('../JSON_input_params/RUN_CNAP_convergence_test_fiber_diameter_precision2.json',...
        parameter_sweep_values,parameter_sweep_variable_name,parameter_sweep_xlabel,title_labels,reference_CNAP_label,x_bounds_unmyelinated);

end

if (flag_run_linear)
    % Linear
    parameter_sweep_xlabel = 'number of templates';
    reference_CNAP_label = '193';
    % reference_signals = arrayfun(@(x) sprintf('../results/ground_truth_CNAP_myel_%02dmm',x),...
    % parameter_sweep_values,'UniformOutput',false);
    run_parameter_sweep_test('../JSON_input_params/RUN_CNAP_convergence_test_number_of_templates_linear1.json',...
        parameter_sweep_values,parameter_sweep_variable_name,parameter_sweep_xlabel,title_labels,reference_CNAP_label,x_bounds_myelinated);
    reference_CNAP_label = '97';
    run_parameter_sweep_test('../JSON_input_params/RUN_CNAP_convergence_test_number_of_templates_linear2.json',...
        parameter_sweep_values,parameter_sweep_variable_name,parameter_sweep_xlabel,title_labels,reference_CNAP_label,x_bounds_unmyelinated);
end
