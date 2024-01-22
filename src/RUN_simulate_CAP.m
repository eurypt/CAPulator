%{

Author: Edgar Pena

This generates the data for and plots every figure in the main manuscript text.
It simulates CAP signals across different conduction distances and fiber
parameters. The only exception is that it does not run the simulations across
tissue conductivity values; instead, the output of those is loaded here and
plotted. To run a desired subset of the simulations/data, edit the
function_run_status_and_handle variable in the main function as indicated by
those comments.

%}
function RUN_simulate_CAP()

% Add the path to the jsonlab repo; this is used to load the JSON files While
% MATLAB has a builtin JSON parser, here we use JSONLab instead of builtin
% jsonencode and jsondecode because JSONLab is an NIH-funded project as of 2021
% that is meant to lead to better standardization across neuroimaging community
addpath('external_dependencies/jsonlab');

% Add the path to the gramm plotting code
addpath('external_dependencies/gramm');

% Specify which analyses to run by setting the number to 0 or 1 next to
% each function handle in the cell array below; a value of 0 will cause
% that function to not be run
function_run_status_and_handle = {
    1, @run_baseline_CNAP_signal % Generate baseline data for Figure 1 & 3
    1, @run_tuned_CNAP_signal % Generatetune data for Figure 11
    1, @plot_comparison_model_vs_in_vivo % Plot Figure 1 & Figure 11
    1, @plot_comparison_brute_force_vs_efficient % Plot Figure 3
    1, @plot_sensitivity_analysis_results % Plot Figures 4 & 5 (as well as some supplementary figures)
    1, @run_conduction_distance_sensitivity_analysis % Generate data for Figure 6
    1, @plot_conduction_distance_sensitivity_analysis % Plot Figure 6
    1, @quantify_random_sampling_effects % Generate data for Figure 7 & plot it
    1, @run_compare_Havton_to_Soltanpour % Generate data for Figure 8
    1, @plot_compare_Havton_to_Soltanpour % Plot Figure 8
    1, @run_evalulate_CV_vs_fiberD_effect_on_CAP % Generate data for Figure 9
    1, @plot_evalulate_CV_vs_fiberD_effect_on_CAP % Plot Figure 9
    1, @quantify_shrinkage_effects % Generate data for Figure 10 & plot it
    };
% Print the date and time this script was run
fprintf('Running script %s at %s\n',mfilename,datetime('now'));
fprintf('MATLAB version: %s\n',version);
for i = 1:length(function_run_status_and_handle)
    run_status = function_run_status_and_handle{i,1};
    function_handle = function_run_status_and_handle{i,2};
    if (run_status)
        % Print which function is running and time it; also print all this info to a file named "diary"
        fprintf('Running function %s\n',func2str(function_handle));
        tic;
        % Run the function
        function_handle();
        % Print the elapsed time
        fprintf('Elapsed time for function %s: %f seconds\n',func2str(function_handle),toc);
        % Skip a line
        fprintf('\n');
    end
end

% Simulate the CAP signals for the 20221004 experiment using the default
% parameters
function run_baseline_CNAP_signal()

% Myelinated fibers
% Specify the parameter files to use
parameter_files = {
    '../JSON_input_params/RUN_CNAP_20221004_ASCENT_myel.json','myelinated',1
    '../JSON_input_params/RUN_CNAP_20221004_ASCENT_unmyel.json','unmyelinated',1
    };

% Save the data to a table. The variables to save are t_window (to be saved as
% time_vector_ms) and ch_neural_avg(:,end) (to be saved as CNAP_uV), as well as
% the channel number, the recording volume conductor, the conduction distance,
% the fiber type, and the data type
model_data_table = [];
for i = 1:size(parameter_files,1)
    params = loadjson(parameter_files{i,1});
    
    % Print a message indicating which parameter file is being used
    fprintf('Simulating  %s\n',parameter_files{i,1});
    % Run the simulation
    [CNAP_signal_uV, final_time_vector_ms,SFAPs_uV] = CAPulator.calculate_CAP_via_temporal_templates(params);
    channel_numbers = parameter_files{i,3};
    conduction_distance = 11 - params.stim_location_change_mm;
    data_type = {'model'};
    fiber_type = {parameter_files{i,2}};
    
    % Save the data to a table
    data_table_i = table({params.extracellular_recording_model_filename},{final_time_vector_ms},{CNAP_signal_uV},...
        channel_numbers,conduction_distance,fiber_type,data_type, ...
        'VariableNames',{'volume_conductor','time_vector_ms','CNAP_uV','channel_number','conduction_distance_mm','fiber_type','data_type'});
    % Append the data table to the in_vivo_data_table
    model_data_table = [model_data_table; data_table_i];
end

% Save model table
save('../results/baseline_CNAP_model_data_table.mat','model_data_table');


% Simulate the CAP signals for the 20221004 experiment using the tuned
% parameters
function run_tuned_CNAP_signal()

% Myelinated fibers
% Specify the parameter files to use
parameter_files = {
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_myel_11mm_ch1_polarity1.json','myelinated',1
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_myel_11mm_ch2_polarity1.json','myelinated',2
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_myel_11mm_ch3_polarity1.json','myelinated',3
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_11mm_ch1_polarity1.json','unmyelinated',1
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_11mm_ch2_polarity1.json','unmyelinated',2
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_11mm_ch3_polarity1.json','unmyelinated',3
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_myel_15mm_ch1_polarity1.json','myelinated',1
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_15mm_ch1_polarity1.json','unmyelinated',1
    '../JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_6mm_ch1_polarity1.json','unmyelinated',1
    };

% Save the data to a table. The variables to save are t_window (to be saved as
% time_vector_ms) and ch_neural_avg(:,end) (to be saved as CNAP_uV), as well as
% the channel number, the recording volume conductor, the conduction distance,
% the fiber type, and the data type
model_data_table = [];
for i = 1:size(parameter_files,1)
    params = loadjson(parameter_files{i,1});
    
    % Print a message indicating which parameter file is being used
    fprintf('Simulating  %s\n',parameter_files{i,1});
    % Run the simulation
    [CNAP_signal_uV, final_time_vector_ms,SFAPs_uV] = CAPulator.calculate_CAP_via_temporal_templates(params);
    channel_numbers = parameter_files{i,3};
    conduction_distance = 11 - params.stim_location_change_mm;
    data_type = {'model'};
    fiber_type = {parameter_files{i,2}};
    
    % Save the data to a table
    data_table_i = table({params.extracellular_recording_model_filename},{final_time_vector_ms},{CNAP_signal_uV},...
        channel_numbers,conduction_distance,fiber_type,data_type, ...
        'VariableNames',{'volume_conductor','time_vector_ms','CNAP_uV','channel_number','conduction_distance_mm','fiber_type','data_type'});
    % Append the data table to the in_vivo_data_table
    model_data_table = [model_data_table; data_table_i];
end

% Save model table
save('../results/tuned_CNAP_model_data_table.mat','model_data_table');


function plot_comparison_model_vs_in_vivo()
% Load the in vivo data
load('../results/in_vivo_data_table.mat','in_vivo_data_table');

% Compare the baseline model data to in vivo; plot all the in vivo data first
% Load the model data to plot on top of the in vivo data
load('../results/baseline_CNAP_model_data_table.mat','model_data_table');
% Plot the results
clearvars g
g = gramm('x',[in_vivo_data_table.time_vector_ms;model_data_table.time_vector_ms],...
    'y',[in_vivo_data_table.CNAP_uV;model_data_table.CNAP_uV],...
    'linestyle',[in_vivo_data_table.data_type;model_data_table.data_type],...
    'subset',[in_vivo_data_table.channel_number;model_data_table.channel_number]==1 &...
    [in_vivo_data_table.conduction_distance_mm;model_data_table.conduction_distance_mm]==11);
g.facet_grid([],[in_vivo_data_table.fiber_type;model_data_table.fiber_type],"scale","independent");
g.geom_line();
g.set_names('x','time (ms)','y','signal (\muV)','column','fiber type');
g.set_text_options('interpreter','tex','base_size',16);
g.set_line_options("styles",{'-',':'});
figure('position',[99,497,803,420]);
g.draw();

% Set the x-axis limits for myelinated fibers
set(g.facet_axes_handles(1),'XLim',[0.16,2.5]);
legend(g.facet_axes_handles(1),{'in vivo','model'})
% Set the x-axis limits for unmyelinated fibers
set(g.facet_axes_handles(2),'XLim',[4,35]);
legend(g.facet_axes_handles(2),{'in vivo','model'})

% Now compare the tuned model data to in vivo; plot all the in vivo data first
clearvars model_data_table
load('../results/tuned_CNAP_model_data_table.mat','model_data_table');
clearvars g
g(1,1) = gramm('x',[in_vivo_data_table.time_vector_ms;model_data_table.time_vector_ms],...
    'y',[in_vivo_data_table.CNAP_uV;model_data_table.CNAP_uV],...
    'color',[in_vivo_data_table.conduction_distance_mm;model_data_table.conduction_distance_mm],...
    'linestyle',[in_vivo_data_table.data_type;model_data_table.data_type],...
    'subset',[in_vivo_data_table.channel_number;model_data_table.channel_number]==1);
g(1,1).facet_grid([in_vivo_data_table.fiber_type;model_data_table.fiber_type],[],"scale","independent",'row_labels',false);
g(1,1).geom_line();
g(1,1).set_names('x','time (ms)','y','signal (\muV)','row','fiber type','color','conduction distance (mm)');
g(1,1).set_text_options('interpreter','tex','base_size',16);
g(1,1).set_layout_options('legend',false);
g(1,1).set_line_options("styles",{'-',':'});
g(1,1).set_color_options('map','matlab');

g(1,2) = gramm('x',[in_vivo_data_table.time_vector_ms;model_data_table.time_vector_ms],...
    'y',[in_vivo_data_table.CNAP_uV;model_data_table.CNAP_uV],...
    'size',[in_vivo_data_table.channel_number;model_data_table.channel_number],...
    'linestyle',[in_vivo_data_table.data_type;model_data_table.data_type],...
    'subset',[in_vivo_data_table.conduction_distance_mm;model_data_table.conduction_distance_mm]==11);
g(1,2).facet_grid([in_vivo_data_table.fiber_type;model_data_table.fiber_type],[],"scale","independent",'row_labels',false);
g(1,2).geom_line();
g(1,2).set_names('x','time (ms)','y','signal (\muV)','row','fiber type');
g(1,2).set_layout_options('legend',false);
g(1,2).set_line_options("styles",{'-',':'});
g(1,2).set_text_options('interpreter','tex','base_size',16);
g(1,2).set_color_options('map',[0.85,0.33,0.10]);

figure('position',[464,240,1101,738]);
g.draw();
set(g(1).facet_axes_handles(1),'XLim',[0.16,2],'YLim',200*[-1 1]);
set(g(1).facet_axes_handles(2),'XLim',[2.5,35],'YLim',3000*[-1 1]);
set(g(2).facet_axes_handles(1),'XLim',[0.16,2],'YLim',200*[-1 1]);
set(g(2).facet_axes_handles(2),'XLim',[2.5,35],'YLim',800*[-1 1]);




function plot_comparison_brute_force_vs_efficient()
% Load the in vivo data
clearvars brute_force
brute_force(1) = load('../results/ground_truth_CNAP_myel.mat','time_vector','CNAP_signal_uV','time_offset_ms');
brute_force(2) = load('../results/ground_truth_CNAP_unmyel.mat','time_vector','CNAP_signal_uV','time_offset_ms');

% Plot the results
clearvars g
time_vectors = cellfun(@(g1,g2) g1-2, {brute_force.time_vector}, {brute_force.time_offset_ms},'UniformOutput',false);
g = gramm('x',time_vectors,...
    'y',{brute_force.CNAP_signal_uV});
g.facet_grid([],{'myelinated';'unmyelinated'},"scale","independent");
g.geom_line();
g.set_names('x','time (ms)','y','signal (\muV)','column','fiber type');
g.set_text_options('interpreter','tex','base_size',16);
g.set_color_options('map','matlab')
figure('position',[593,211,1125,471]);
g.draw();

% Load the model data to plot on top of the in vivo data
load('../results/baseline_CNAP_model_data_table.mat','model_data_table');
g.update('x',model_data_table.time_vector_ms,...
    'y',model_data_table.CNAP_uV,...
    'subset',model_data_table.channel_number==1 & model_data_table.conduction_distance_mm==11);
g.geom_line();
g.set_line_options("styles",{':'});
g.set_color_options('map',[0,0,0])
g.draw();

% Set the x-axis limits for myelinated fibers
set(g.facet_axes_handles(1),'XLim',[0.16,2.5]);
legend(g.facet_axes_handles(1),{'brute force','filtered interpolated templates'})
% Set the x-axis limits for unmyelinated fibers
set(g.facet_axes_handles(2),'XLim',[4,35]);
legend(g.facet_axes_handles(2),{'brute force','filtered interpolated templates'})


function run_conduction_distance_sensitivity_analysis()

for i = 1:2
    if (i==1)
        base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_myel.json';
        workspace_filename = '../results/conduction_distance_effects_myel';
    elseif (i==2)
        base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_unmyel.json';
        workspace_filename = '../results/conduction_distance_effects_unmyel';
    else
        error('invalid value of i: %d')
    end
    
    % Load base parameters file
    params = loadjson(base_parameters_filename);
    
    % Run the CNAP after changing the multiple conduction distances
    n_distances = 10;
    % the default edge-to-edge distance is 9232.5 um (when the center-to-center distance is 11 mm with a
    % bipolar microleads stim cuff and a tripolar microleads rec cuff).
    % Adjust this distance to be 0 mm edge-to-edge (-9232.5 um) and then
    % add whatever the desired conduction distance is
    % So add 9232.5 um, which corresponds to an edge-to-edge distance is 0 mm
    %     stim_location_change_mm = 9232.5 - logspace(log10(1e3),log10(100e3),n_distances);
    KNOWN_CONDUCTION_DISTANCE_MM = 11; % [mm]
    %
    edge_to_edge_distance_mm = logspace(log10(4),log10(100),n_distances);
    stim_location_change_mm = KNOWN_CONDUCTION_DISTANCE_MM - 1.7675 - edge_to_edge_distance_mm;
    
    % indicate the known location of the recording electrode to later label
    % the plots in terms of conduction distance rather than AP start
    % location
    output_table = [];
    for idx = 1:length(stim_location_change_mm)
        % Print a message indicating which parameter file is being used
        fprintf('Simulating conduction distance %d of %d with %s\n',idx,length(stim_location_change_mm),base_parameters_filename);
        
        % set first AP location
        params.stim_location_change_mm = stim_location_change_mm(idx);
        output_table = vertcat(output_table,process_params(params));
    end
    
    % calculate conduction distance to label the plots
    conduction_distances_all_mm = KNOWN_CONDUCTION_DISTANCE_MM - stim_location_change_mm;
    
    %{
        The workspace variables include the
        following key variables:
        - output_table: table of the CNAP signals for each conduction distance.
            This table contains the following columns: common_time_vector_ms,
            CAP_signal_uV, params
        - conduction_distances_all_mm: conduction distance for each signal
        - common_time_vector_ms: common time vector for each signal
        - params (the input parameters)
        - edge_to_edge_distance_mm: edge-to-edge distance for each signal
        - stim_location_change_mm: stim location change for each signal
        - KNOWN_CONDUCTION_DISTANCE_MM: known conduction distance
    %}
    save(workspace_filename);
end

function plot_conduction_distance_sensitivity_analysis()

output_table = [];
for i = 1:2
    if (i==1)
        workspace_filename = '../results/conduction_distance_effects_myel';
        fiber_type = 'myelinated';
    elseif (i==2)
        workspace_filename = '../results/conduction_distance_effects_unmyel';
        fiber_type = 'unmyelinated';
    else
        error('invalid value of i: %d')
    end
    
    dataset_i = load(workspace_filename,'output_table','conduction_distances_all_mm');
    % Add a column to the table to indicate the fiber type
    dataset_i.output_table.fiber_type = repmat({fiber_type},size(dataset_i.output_table,1),1);
    % Add a column to the table to indicate the conduction distance; conduction distance has the same number of rows as the table
    dataset_i.output_table.conduction_distances_all_mm = dataset_i.conduction_distances_all_mm';
    output_table = vertcat(output_table,dataset_i.output_table);
    
end

% Plot subset of signals to illustrate temporal dispersion
clearvars g
% only a subset of conduction_distances_all_mm values for visual clarity
g = gramm('x',output_table.common_time_vector_ms,...
    'y',cellfun(@(signal) signal/1e3, output_table.CAP_signal_uV, 'UniformOutput',false),...
    'color',round(output_table.conduction_distances_all_mm,1),...
    'subset',ismember(round(output_table.conduction_distances_all_mm,1),[5.8,9.9,18.5]));
g.geom_line();
g.fig(output_table.fiber_type);
g.set_title('');
g.axe_property('XLim',[0 60]);
%     g.set_continuous_color('colormap','parula')
g.set_names('x','time (ms)','y','signal (mV)','color','conduction distance (mm)');
g.set_text_options('interpreter','tex','base_size',16,...
    'legend_scaling',0.8,'legend_title_scaling',0.9);
g.set_layout_options('legend_position',[0.55 0.65 0.45 0.35]);
figure;
g.draw();
% Set the x-axis limits for myelinated fibers
set(g.facet_axes_handles(1),'XLim',[0,2.5]);

% Plot summary of amplitude for all of the signals
clearvars g
g = gramm('x',output_table.conduction_distances_all_mm,...
    'y',cellfun(@(signal) max(signal/1e3)-min(signal/1e3), ...
    output_table.CAP_signal_uV),'color',output_table.fiber_type);
g.geom_line();
g.geom_point();
g.axe_property('XScale','log','YScale','log',...
    'XLim',[min(output_table.conduction_distances_all_mm),max(output_table.conduction_distances_all_mm)],...
    'YLim',[0.08 22]);
g.set_names('x','conduction distance (mm)','y','V_{pk-pk} (mV)','color','fiber type');
g.set_text_options('interpreter','tex','base_size',16,...
    'legend_scaling',0.8,'legend_title_scaling',0.9);
g.set_layout_options('legend_position',[0.55 0.6 0.45 0.35]);
g.set_color_options('map','matlab');
figure;
g.draw();



% Plot the nonlinear fit of  the data point such that the log of the fit is a line
clearvars g
g = gramm('x',(output_table.conduction_distances_all_mm),...
    'y',cellfun(@(signal) (max(signal/1e3)-min(signal/1e3)), ...
    output_table.CAP_signal_uV),'subset',output_table.conduction_distances_all_mm<26,'color',output_table.fiber_type);
g.stat_fit('fun',@(param1,param2,x)param2*x.^param1,'intopt','functional');
g.axe_property('XScale','log','YScale','log',...
    'XLim',[min(output_table.conduction_distances_all_mm),max(output_table.conduction_distances_all_mm)],...
    'YLim',[0.08 22]);
g.set_names('x','conduction distance (mm)','y','V_{pk-pk} (mV)');
g.set_text_options('interpreter','tex','base_size',16,...
    'legend_scaling',0.8,'legend_title_scaling',0.9);
g.set_layout_options('legend_position',[0.55 0.6 0.45 0.35]);
% set line style to be dotted
g.set_line_options('styles',{'--'});
figure;
g.draw();
% Print the stats info for each fiber type
g.results.stat_fit.model


% Compare CNAP from the Havton fiber diameters to the CNAP from the
% Soltanpour distribution
function run_compare_Havton_to_Soltanpour()

% run comparisons for both myelinated and ummyelinated
output_data_structure.flag_myelinated = [];
output_data_structure.pseudonym = {};
output_data_structure.CAP_signal_uV = {};
output_data_structure.common_time_vector_ms = {};
output_data_structure.fiber_diameters_in_nerve_um = {};

% Edit the fiber diameters to be cVN & run CNAP simulation
for i = 1:4
    switch i
        case 1
            flag_myelinated = 1;
            base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_myel.json';
            pseudonym = 'cVN';
        case 2
            flag_myelinated = 0;
            base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_unmyel.json';
            pseudonym = 'cVN';
        case 3
            flag_myelinated = 1;
            base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_myel_Soltanpour.json';
            pseudonym = 'cVN_Soltanpour';
        case 4
            flag_myelinated = 0;
            base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_unmyel_Soltanpour.json';
            pseudonym = 'cVN_Soltanpour';
        otherwise
            error('invalid value of i: %d')
    end
    
    % Load base parameters file
    params = loadjson(base_parameters_filename);
    
    [CAP_signal_uV, common_time_vector_ms, ~] = CAPulator.calculate_CAP_via_temporal_templates(params);
    
    output_data_structure.pseudonym = vertcat(...
        output_data_structure.pseudonym, {pseudonym});
    output_data_structure.flag_myelinated = vertcat(...
        output_data_structure.flag_myelinated, flag_myelinated);
    output_data_structure.CAP_signal_uV = vertcat(...
        output_data_structure.CAP_signal_uV,{CAP_signal_uV});
    output_data_structure.common_time_vector_ms = vertcat(...
        output_data_structure.common_time_vector_ms,{common_time_vector_ms});
    output_data_structure.fiber_diameters_in_nerve_um = vertcat(...
        output_data_structure.fiber_diameters_in_nerve_um,{params.fiber_diameters_in_nerve_um});
end

% Save the workspace
save('../results/compare_Havton_to_Soltanpour_workspace.mat');

function plot_compare_Havton_to_Soltanpour()

% Load the workspace
load('../results/compare_Havton_to_Soltanpour_workspace.mat');

for i = 0:1
    % Set edges based on Soltanpour & Santer 1996 data such that the bin centers
    % from those data are the bin centers for the specified edges
    if (i==0) % unmyelinated
        x_bounds = [0,24]; % [ms]
        edges = 0.21:0.14:1.33; % [um]
    else % myelinated
        x_bounds = [0,2.5]; % [ms]
        edges = 0.5:1:8.5; % [um]
    end
    % Plot CNAPs
    clearvars g
    g = gramm('x',output_data_structure.common_time_vector_ms,...
        'y',cellfun(@(signal) signal/1e3, output_data_structure.CAP_signal_uV, 'UniformOutput',false),...
        'color',output_data_structure.pseudonym,'subset',output_data_structure.flag_myelinated==i);
    % disable column labels
    g.set_order_options('column',-1);
    g.facet_grid([],output_data_structure.flag_myelinated,'scale','independent','column_labels',false);
    % set the x-axis based on the column facet value
    g.axe_property('XLim',x_bounds);
    g.set_names('x','time (ms)','y','signal (mV)','color','source');
    g.set_text_options('interpreter','tex','base_size',16,...
        'legend_scaling',0.8,'legend_title_scaling',0.9);
    g.geom_line();
    figure('position',[680,679,885,299]);
    g.draw();
    
    % Plot the histograms of fiber diameters
    clearvars g
    g = gramm('x',output_data_structure.fiber_diameters_in_nerve_um,...
        'color',output_data_structure.pseudonym,'subset',output_data_structure.flag_myelinated==i);
    g.set_names('x','fiber diameter (\mum)','color','source');
    g.set_text_options('interpreter','tex','base_size',16,...
        'legend_scaling',0.8,'legend_title_scaling',0.9);
    g.stat_bin('edges',edges,'geom','overlaid_bar');
    figure('position',[680,679,885,299]);
    g.draw();
    
    
end

% Print out summary metrics for each fiber type and the fiber type and pseudonym
% for each median
for ind = 1:length(output_data_structure.fiber_diameters_in_nerve_um)
    % Median
    fprintf('median fiber diameter for %s (myelinated=%d): %0.2f um\n',...
        output_data_structure.pseudonym{ind},...
        output_data_structure.flag_myelinated(ind),...
        median(output_data_structure.fiber_diameters_in_nerve_um{ind}));
    % Total number of fibers
    fprintf('total number of fibers for %s (myelinated=%d): %d\n',...
        output_data_structure.pseudonym{ind},...
        output_data_structure.flag_myelinated(ind),...
        length(output_data_structure.fiber_diameters_in_nerve_um{ind}));
end


%{

Calculate CNAPs using  different fiberD-to-CV relationships iterate through the
effect of different possible relationships.

%}
function run_evalulate_CV_vs_fiberD_effect_on_CAP()
random_seeds = [10 15 33];

% define the possible options for jitter in the z direction; 0 means no
% jitter; 1 means jitter
z_jitter_options = [0 1];

% list the following:
% model name (for labeling), conduction distance, myel flag, random seed, and function of fiberD vs CV

all_pseudonyms_and_distances = {
    % same as model A001, except D vs. V is fit to the biophysical model D vs. V
    'model default',14e3,1,random_seeds(1),@(D,V) 3.877*D - 2.491 ,z_jitter_options(1),1
    
    % same as model A002, except slope coefficient is halved/double
    'faster',14e3,1,random_seeds(1),@(D,V) (3.877*sqrt(2))*D - 2.491 ,z_jitter_options(1),1
    'slower',14e3,1,random_seeds(1),@(D,V) (3.877/sqrt(2))*D - 2.491 ,z_jitter_options(1),1
    
    % same as model A002, except D vs. V is fit to Hursh 1939
    'Hursh 1939',14e3,1,random_seeds(1),@(D,V) 5.972*D - 3.297,z_jitter_options(1),1
    
    % same as model A007, except D vs. V is fit to the biophysical model D vs. V
    'model default',14e3,0,random_seeds(1),@(D,V) 0.7842*sqrt(D) + 0.004242,z_jitter_options(1),1
    
    % same as model A008, except slope coefficient is halved/double
    'faster',14e3,0,random_seeds(1),@(D,V) (0.7842*sqrt(2))*sqrt(D) + 0.004242,z_jitter_options(1),1
    'slower',14e3,0,random_seeds(1),@(D,V) (0.7842/sqrt(2))*sqrt(D) + 0.004242,z_jitter_options(1),1
    
    % same as model A008, except D vs. V is linear fit to Hoffmeister 1991 Fig 5A
    'Hoffmeister 1991',14e3,0,random_seeds(1),@(D,V) 0.6496*D + 0.1309,z_jitter_options(1),1
    
    };


for params_i = 1:size(all_pseudonyms_and_distances,1)
    % Specify pseudonym
    pseudonym = all_pseudonyms_and_distances{params_i,1};
    
    flag_myelinated = all_pseudonyms_and_distances{params_i,3};
    
    % manually modify the fiberD-to-CV relationship in a temporary template
    % data source file to iterate through the effect of different possible
    % relationships.
    if (flag_myelinated)
        base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_myel.json';
        fiber_type_label = 'myel';
        adjust_ms = 0.1; % [ms] adjust for the stimulus artifact
    else
        base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_unmyel.json';
        fiber_type_label = 'unmyel';
        adjust_ms = 0.4; % [ms] adjust for the stimulus artifact
    end
    
    % Load base parameters file
    params = loadjson(base_parameters_filename);
    
    % Load the original template data
    original_template_data = load(params.template_data_source_filename,'output_data_structure','fiber_type');
    
    % Set new CV values for each fiber diameter using a new specified
    % relationship
    old_CV = vertcat([original_template_data.output_data_structure.all_conduction_velocity_m_per_s]);
    new_CV = all_pseudonyms_and_distances{params_i,5}(...
        vertcat([original_template_data.output_data_structure.all_fiber_diameters]),...
        old_CV);
    % Also adjust reference_peak_times_all_templates to reflect the new CV
    assert(all(new_CV>0),'All new CV values must be positive')
    for cv_ind = 1:length(original_template_data.output_data_structure)
        original_template_data.output_data_structure(cv_ind).reference_peak_time_i_dipole_ms = ...
            (original_template_data.output_data_structure(cv_ind).reference_peak_time_i_dipole_ms+adjust_ms) * ...
            old_CV(cv_ind) / new_CV(cv_ind) - adjust_ms;
        original_template_data.output_data_structure(cv_ind).reference_peak_time_i_ms = ...
            (original_template_data.output_data_structure(cv_ind).reference_peak_time_i_ms+adjust_ms) * ...
            old_CV(cv_ind) / new_CV(cv_ind) - adjust_ms;
        original_template_data.output_data_structure(cv_ind).all_conduction_velocity_m_per_s = new_CV(cv_ind);
        
    end
    
    % Specify the temporary template data filename; the modified CV values and
    % original other data will be stored in here for CAP generation
    temporary_template_data_source_filename = 'temporary_template_data.mat';
    save(temporary_template_data_source_filename,'-struct','original_template_data');
    
    % Update the params structure to point to the temporary template data
    params.template_data_source_filename = temporary_template_data_source_filename;
    
    % Calculate CAP using temporal template method
    [CAP_signal_uV, common_time_vector_ms, ~] = CAPulator.calculate_CAP_via_temporal_templates(params);
    
    % Store CNAP data into a structure
    if (params_i==1)
        output_data_structure = [];
        output_data_structure.pseudonym = {pseudonym};
        output_data_structure.flag_myelinated = flag_myelinated;
        output_data_structure.CAP_signal_uV = {CAP_signal_uV};
        output_data_structure.common_time_vector_ms = {common_time_vector_ms};
    else
        output_data_structure.pseudonym = vertcat(...
            output_data_structure.pseudonym, {pseudonym});
        output_data_structure.flag_myelinated = vertcat(...
            output_data_structure.flag_myelinated,flag_myelinated);
        output_data_structure.CAP_signal_uV = vertcat(...
            output_data_structure.CAP_signal_uV,{CAP_signal_uV});
        output_data_structure.common_time_vector_ms = vertcat(...
            output_data_structure.common_time_vector_ms,{common_time_vector_ms});
    end
    
    % Clean up: Delete the temporary template data file
    delete(temporary_template_data_source_filename);
    
end

% Save output_data_structure
save('../results/CV_vs_fiberD_effect.mat','-struct','output_data_structure');

function plot_evalulate_CV_vs_fiberD_effect_on_CAP()
output_data_structure = load('../results/CV_vs_fiberD_effect.mat');

% Plot the CAP signals
for i = 0:1
    if (i==0) % unmyelinated
        x_bounds = [0,40]; % [ms]
        svg_filename = 'CV_vs_D_effect_myelinated';
    else % myelinated
        x_bounds = [0,2.5]; % [ms]
        svg_filename = 'CV_vs_D_effect_unmyelinated';
    end
    
    clearvars g
    g = gramm('x',output_data_structure.common_time_vector_ms,...
        'y',cellfun(@(signal) signal/1e3, output_data_structure.CAP_signal_uV, 'UniformOutput',false),...
        'color',output_data_structure.pseudonym,...
        'subset',output_data_structure.flag_myelinated==i);
    g.geom_line();
    g.set_names('x','Time (ms)','y','CAP (mV)','color','Model');
    g.set_text_options('interpreter','tex','base_size',16,...
        'legend_scaling',0.8,'legend_title_scaling',0.9);
    g.axe_property('XLim',x_bounds);
    figure;
    g.draw();
    % g.export('file_name',svg_filename,'file_type','svg');
end


% Quantify the effects of fiber shrinkage
function quantify_shrinkage_effects()

% run comparisons for both myelinated and ummyelinated
for i = 1:2
    if (i==1)
        flag_myelinated = 1;
        base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_myel.json';
        x_bounds = [0 2];
        time_series_fig_filename = '../results/shrinkage_correction_effect_myelinated';
        summary_fig_filename = '../results/shrinkage_correction_effect_myelinated_summary';
    elseif (i==2)
        flag_myelinated = 0;
        base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_unmyel.json';
        x_bounds = [0 25];
        time_series_fig_filename = '../results/shrinkage_correction_effect_unmyelinated';
        summary_fig_filename = '../results/shrinkage_correction_effect_unmyelinated_summary';
    else
        error('invalid value of i: %d')
    end
    
    % Load base parameters file
    params = loadjson(base_parameters_filename);
    
    % Run the CNAP across multiple shrinkage correction factors
    shrinkage_correction_factors = [1:0.05:1.2];
    original_fiber_diameters_um = params.fiber_diameters_in_nerve_um;
    output_table = [];
    for idx=1:length(shrinkage_correction_factors)
        params.fiber_diameters_in_nerve_um = shrinkage_correction_factors(idx)*original_fiber_diameters_um;
        output_table = vertcat(output_table,process_params(params));
    end
    
    % Plot all CNAPs vs. time (overlaid)
    figure('color',[1 1 1]);
    for idx=1:length(shrinkage_correction_factors)
        plot(output_table.common_time_vector_ms{idx},1e-3*output_table.CAP_signal_uV{idx});
        xlabel('time (ms)');
        ylabel('CNAP (mV)');
        set(gca,'FontSize',14);
        hold on;
    end
    xlim(x_bounds);
    % change the line colors so that they are unique and more visible
    lines = get(gca,'Children');
    new_colors = parula(length(lines));
    for line_ind = 1:length(lines)
        lines(line_ind).Color = new_colors(line_ind,:);
    end
    legend(arrayfun(@(x) num2str(x,'%0.02f'), shrinkage_correction_factors,'UniformOutput',false),...
        'location','best','numcolumns',2);
    savefig(gcf,time_series_fig_filename);
    
    % Calculate and plot the peak-to-peak voltage and the latency
    Vpk2pk_mV = [];
    latency = [];
    for idx = 1:length(shrinkage_correction_factors)
        Vpk2pk_mV(idx) = 1e-3*(max(output_table.CAP_signal_uV{idx})-min(output_table.CAP_signal_uV{idx}));
        [~,negative_peak_ind] = min(output_table.CAP_signal_uV{idx});
        latency(idx) = output_table.common_time_vector_ms{idx}(negative_peak_ind);
    end
    
    % Plot summary figure
    figure('color',[1 1 1]);
    % amplitude summary
    yyaxis left
    plot(shrinkage_correction_factors,Vpk2pk_mV,'.-','MarkerSize',20);
    xlabel('shrinkage correction factor');
    ylabel('V_{pk-pk} (mV)');
    set(gca,'FontSize',14);
    % latency summary
    yyaxis right
    plot(shrinkage_correction_factors,latency,'.-','MarkerSize',20);
    xlabel('shrinkage correction factor');
    ylabel('negative peak latency (ms)');
    set(gca,'FontSize',14);
    
    % save figure
    savefig(gcf,summary_fig_filename);
end

function quantify_random_sampling_effects()


% run comparisons for both myelinated and ummyelinated
ouptut_data_structure = [];
for i=1:2
    
    % Load base parameters file
    if (i==1)
        flag_myelinated = 1;
        base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_myel.json';
        x_bounds = [0 3.7];
        y_bounds = 950*[-1 1];
        % I chose the bin sizes to include 0.000001 um
        % (ground truth) and to also include a progression from very coarse
        % (~2 um) with roughly doubling precision until a precision that
        % was basically convered with the ground truth
        n_bins = [8796306,4,9,18,36,64,128];
        random_samples_filename = 'random_samples_myelinated';
        nbins_effects_filename = 'nbins_effects_myelinated';
        histograms_nbins_filename = 'histograms_nbins_myelinated';
        
    elseif (i==2)
        flag_myelinated = 0;
        base_parameters_filename = '../JSON_input_params/RUN_CNAP_20221004_ASCENT_unmyel.json';
        x_bounds = [0 35];
        y_bounds = 4700*[-1 1];
        % I chose the bin sizes to include 0.000001 um
        % (ground truth) and to also include a progression from very coarse
        % (~0.28 um) with roughly doubling precision until a precision that
        % was basically convered with the ground truth
        % Soltanpour & Santer used 0.14 um for their histogram, so the goal
        % was to include approximately that precision in the range
        n_bins = [1790945,6,13,26,52,104];
        random_samples_filename = 'random_samples_unmyelinated';
        nbins_effects_filename = 'nbins_effects_unmyelinated';
        histograms_nbins_filename = 'histograms_nbins_unmyelinated';
        
        
    else
        error('invalid value of i: %d')
    end
    
    params = loadjson(base_parameters_filename);
    
    % Define base fiber diameters
    original_values = params.fiber_diameters_in_nerve_um;
    
    
    % Run the CNAP across multiple random sets of fiber diameters derived
    % form the same fiber diameter distribution
    bin_derived_values = {};
    counts = {};
    centers = {};
    bin_sizes = [];
    for bin_ind=1:length(n_bins)
        
        edges_i = linspace(min(original_values),max(original_values),n_bins(bin_ind)+1);
        [counts_i,~,bin] = histcounts(original_values,edges_i);
        bin_sizes(bin_ind) = mode(diff(edges_i));
        centers_i = edges_i(1:(end-1)) + bin_sizes(bin_ind)/2;
        bin_derived_values{bin_ind} = centers_i(bin);
        % add one more edge above and below the chosen edges so that the
        % zero counts are shown in those during plotting
        edges_i = [edges_i(1)-bin_sizes(bin_ind),edges_i,edges_i(end)+bin_sizes(bin_ind)];
        counts{bin_ind} = [0,counts_i,0];
        centers{bin_ind} = [min(centers_i)-bin_sizes(bin_ind),centers_i,max(centers_i)+bin_sizes(bin_ind)];
    end
    
    % make a pseudo-histogram plot using the centers and counts by making a
    % vector of edges; I did it this way because the gramm overlaid_bar
    % option does not allow changing the bin size for each group
    edges = cellfun(@(c,s) reshape([c-s/2; c+s/2],1,[]), centers,num2cell(bin_sizes),'UniformOutput',false);
    edge_counts = cellfun(@(t) reshape([t; t],1,[]), counts,'UniformOutput',false);
    clearvars g
    g = gramm('x',edges,'y',edge_counts,'color',round(bin_sizes,2));
    g.geom_line();
    %     g.stat_bin('geom','overlaid_bar','width',bin_sizes);
    g.set_names('x','fiber diameter (\mum)','y','# fibers','color','bin size (\mum)');
    g.set_text_options('interpreter','tex','base_size',14,'title_scaling',1);
    g.facet_grid([],bin_sizes,'scale','independent','column_labels',false);
    g.set_layout_options('legend',true);
    g.axe_property('XLim',[0,max(original_values)]); %,'YLim',700*[-1 1]);
    figure('position',[244,697,1526,246]);
    g.draw();
    % g.export('file_name',histograms_nbins_filename,'file_type','svg');
    
    
    output_table = [];
    for bin_ind=1:length(n_bins)
        params.fiber_diameters_in_nerve_um = bin_derived_values{bin_ind}';
        output_table_i = process_params(params);
        output_table = vertcat(output_table,output_table_i);
    end
    
    
    clearvars g
    g = gramm('x',output_table.common_time_vector_ms,'y',output_table.CAP_signal_uV,...
        'color',round(bin_sizes,2));
    g.geom_line();
    g.set_names('x','time (ms)','y','signal (\muV)','color','bin size (\mum)',...
        'column','bin size (\mum)');
    g.set_text_options('interpreter','tex','base_size',14,'title_scaling',1);
    g.facet_grid([],bin_sizes,'column_labels',false);
    g.set_layout_options('legend',true);
    g.axe_property('XLim',x_bounds,'YLim',y_bounds);
    figure('position',[244,697,1526,281]);
    g.draw();
    % g.export('file_name',nbins_effects_filename,'file_type','svg');
    
    % run each bin size with multiple random seeds
    output_data_structure.random_seed = [];
    output_data_structure.bin_size = [];
    output_data_structure.CAP_signal_uV = {};
    output_data_structure.common_time_vector_ms = {};
    random_seeds = [1 2 3];
    for bin_ind = 1:length(n_bins)
        for seed_ind = 1:length(random_seeds)
            output_data_structure.random_seed = vertcat(output_data_structure.random_seed,...
                random_seeds(seed_ind));
            output_data_structure.bin_size = vertcat(output_data_structure.bin_size,...
                bin_sizes(bin_ind));
            
            %%% Calculate the random samples using inverse transform sampling
            randomly_sampled_fiber_diameters = inverse_transform_sample_histogram(...
                centers{bin_ind},counts{bin_ind},random_seeds(seed_ind),...
                length(bin_derived_values{bin_ind}));
            
            params.fiber_diameters_in_nerve_um = randomly_sampled_fiber_diameters';
            
            [CAP_signal_uV, common_time_vector_ms, ~] = ...
                CAPulator.calculate_CAP_via_temporal_templates(params);
            
            output_data_structure.CAP_signal_uV = vertcat(...
                output_data_structure.CAP_signal_uV,{CAP_signal_uV});
            output_data_structure.common_time_vector_ms = vertcat(...
                output_data_structure.common_time_vector_ms,{common_time_vector_ms});
        end
    end
    
    clearvars g
    %     g = gramm('x',output_data_structure.common_time_vector_ms,'y',output_data_structure.CAP_signal_uV,...
    %         'color',round(output_data_structure.bin_size,2),'subset',output_data_structure.random_seed==1);
    g = gramm('x',output_data_structure.common_time_vector_ms,'y',output_data_structure.CAP_signal_uV,...
        'color',round(output_data_structure.bin_size,2),'linestyle',output_data_structure.random_seed);
    g.geom_line();
    g.set_names('x','time (ms)','y','signal (\muV)','linestyle','random seed',...
        'color','bin size (\mum)');
    g.set_text_options('interpreter','tex','base_size',14,'title_scaling',1);
    g.facet_grid([],output_data_structure.bin_size,'column_labels',false);
    g.set_layout_options('legend',true);
    g.axe_property('XLim',x_bounds,'YLim',y_bounds);
    figure('position',[244,697,1526,281]);
    g.draw();
    % g.export('file_name',random_samples_filename,'file_type','svg');
end

%{

Plot the sensitivity analysis results for the myelinated and unmyelinated

%}
function plot_sensitivity_analysis_results

flag_myelinated = 1;
plot_for_fiber_type(flag_myelinated)

flag_myelinated=0;
plot_for_fiber_type(flag_myelinated)

function plot_for_fiber_type(flag_myelinated)

arguments
    flag_myelinated
end

output_table = load_sensitivity_analysis_results(flag_myelinated);
sensitivity_analysis_table = convert_COMSOL_sweep_txt_to_table('../results/sens_analysis_full.txt');

% Create a data structure that contains the same fields as the table, but
% where each field is parsed for a number, then add in a field containing
% the desired data of interest
sensitivity_analysis_struct = [];
for j = 1:length(sensitivity_analysis_table.Properties.VariableNames)
    for i = 1:length(sensitivity_analysis_table.(sensitivity_analysis_table.Properties.VariableNames{j}))
        sensitivity_analysis_struct.(sensitivity_analysis_table.Properties.VariableNames{j})(i) = ...
            sensitivity_analysis_table.(sensitivity_analysis_table.Properties.VariableNames{j})(i);
    end
end
for i = 1:size(sensitivity_analysis_table,1)
    sensitivity_analysis_struct.anisotropy_ratio_vals(i) = ...
        round(sensitivity_analysis_struct.sigma_endoneurium_z(i)./...
        sensitivity_analysis_struct.sigma_endoneurium_x(i),2);
    
    % Add the desired metric to plot
    sensitivity_analysis_struct.Vpk2pk(i) = max(output_table.CAP_signal{i}) - ...
        min(output_table.CAP_signal{i});
    
    
end
sensitivity_analysis_struct.CAP_signals_all = cell2mat(output_table.CAP_signal');
sensitivity_analysis_struct.common_time_vector_ms = output_table.common_time_vector_ms(1,:)';

%% Plot just a subset of the data to highlight key pieces of the story
if (flag_myelinated)
    ytick_locations = 2.^[5:12];
else
    ytick_locations = 2.^[7:14];
end
for i = 1:4
    switch i
        case 1
            % Amplitude vs. fill
            x_data = sensitivity_analysis_struct.sigma_fill;
            x_data_label = '\sigma_{surround} (S/m)';
            x_tick_label_format = '%0.2f';
            subset_indices = round(sensitivity_analysis_struct.anisotropy_ratio_vals,2)==3.43 & ...
                round(sensitivity_analysis_struct.sigma_endoneurium_z,5)==0.57143 & ...
                round(sensitivity_analysis_struct.sigma_perineurium,6)==round(8.70322e-4,6);
        case 2
            % Amplitude vs. perineurial conductivity
            x_data = sensitivity_analysis_struct.sigma_perineurium;
            x_data_label = '\sigma_{perineurium} (S/m)';
            x_tick_label_format = '%0.1e';
            subset_indices = round(sensitivity_analysis_struct.anisotropy_ratio_vals,2)==3.43 & ...
                round(sensitivity_analysis_struct.sigma_endoneurium_z,5)==0.57143 & ...
                round(sensitivity_analysis_struct.sigma_fill,6)==round(1.76,6);
        case 3
            % Amplitude vs. longitudinal endoneural conductivity
            x_data = sensitivity_analysis_struct.sigma_endoneurium_z;
            x_data_label = '\sigma_{z} (S/m)';
            x_tick_label_format = '%0.2f';
            subset_indices = round(sensitivity_analysis_struct.anisotropy_ratio_vals,2)==3.43 & ...
                round(sensitivity_analysis_struct.sigma_fill,2)==1.76 & ...
                round(sensitivity_analysis_struct.sigma_perineurium,6)==round(8.70322e-4,6);
        case 4
            % Amplitude vs. endoneural anisotropy
            x_data = sensitivity_analysis_struct.anisotropy_ratio_vals;
            x_data_label = '\sigma_z / \sigma_r';
            x_tick_label_format = '%0.2f';
            subset_indices = round(sensitivity_analysis_struct.sigma_endoneurium_z,5)==0.57143 & ...
                round(sensitivity_analysis_struct.sigma_fill,2)==1.76 & ...
                round(sensitivity_analysis_struct.sigma_perineurium,6)==round(8.70322e-4,6);
        otherwise
            error('Invalid case in switch statement');
    end
    
    y_data = 1e6*sensitivity_analysis_struct.Vpk2pk;
    color_data = 360-sensitivity_analysis_struct.Theta_CylUm300t_20230323_001;
    
    clear g
    g=gramm('x',x_data,...
        'y',y_data,...
        'color',color_data,...
        'subset',subset_indices);
    
    g.geom_line();
    g.geom_point();
    
    xtick_vals =  unique(x_data);
    xtick_labels = arrayfun(@(x) num2str(x,x_tick_label_format),xtick_vals,'UniformOutput',false);
    g.set_names('x',x_data_label,'y','V_{pk-pk} (\muV)','color','cuff opening (\deg)');
    g.axe_property('YScale','log','XScale','log','YMinorTick','off','YGrid','on',...
        'YMinorGrid','off','YTick',ytick_locations,'YLim',[min(ytick_locations),max(ytick_locations)],...
        'XMinorTick','off','XLim',...
        [min(x_data),max(x_data)],'XTick',xtick_vals,'XTickLabel',xtick_labels);
    g.set_layout_options('legend',false);
    g.set_text_options('interpreter','tex','base_size',12);
    figure;
    g.draw();
end

%% Plot all the data as a matrix
if (flag_myelinated)
    ytick_locations = 4.^[1:7];
else
    ytick_locations = 4.^[2:8];
end
for i =1:2
    switch i
        case 1
            x_data = sensitivity_analysis_struct.sigma_fill;
            x_data_label = '\sigma_{surround} (S/m)';
            x_tick_label_format = '%0.2f';
            subset_indices = sensitivity_analysis_struct.Theta_contact_CylUm300t_20230323_001==360;
            
        case 2
            x_data = sensitivity_analysis_struct.sigma_fill;
            x_data_label = '\sigma_{surround} (S/m)';
            x_tick_label_format = '%0.2f';
            subset_indices = sensitivity_analysis_struct.Theta_contact_CylUm300t_20230323_001==344;
        otherwise
            error('Invalid case value')
    end
    
    y_data = 1e6*sensitivity_analysis_struct.Vpk2pk;
    color_data = sensitivity_analysis_struct.sigma_perineurium;
    
    clear g
    g=gramm('x',x_data,...
        'y',y_data,...
        'color',color_data,...
        'subset',subset_indices);
    
    g.facet_grid(round(sensitivity_analysis_struct.anisotropy_ratio_vals,2),...
        round(sensitivity_analysis_struct.sigma_endoneurium_z,2));
    g.geom_line();
    g.geom_point();
    
    xtick_vals =  unique(x_data);
    xtick_labels = arrayfun(@(x) num2str(x,x_tick_label_format),xtick_vals,'UniformOutput',false);
    g.set_names('row','\sigma_z / \sigma_r','column','\sigma_z (S/m)','x',x_data_label,...
        'y','V_{pk-pk} (\muV)','color','\sigma_{perineurium} (S/m)');
    g.axe_property('YScale','log','XScale','log','YMinorTick','off','YGrid','on',...
        'YMinorGrid','off','YTick',ytick_locations,'YLim',[min(ytick_locations),max(ytick_locations)],...
        'XMinorTick','off','XLim',...
        [min(x_data),max(x_data)],'XTick',xtick_vals,'XTickLabel',xtick_labels);
    g.set_layout_options('legend',true);
    g.set_text_options('interpreter','tex','base_size',10);
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    g.draw();
end


%% Plot waveforms for a subset of the data to highlight key pieces of the waveform effects story
if (flag_myelinated)
    x_bounds = [0,2];
else
    x_bounds = [2,37];
end
for i = 1:4
    switch i
        case 1
            % Waveform vs. fill
            x_data = sensitivity_analysis_struct.common_time_vector_ms';
            x_data_label = 'time (ms)';
            subset_indices = round(sensitivity_analysis_struct.anisotropy_ratio_vals,2)==3.43 & ...
                round(sensitivity_analysis_struct.sigma_endoneurium_z,5)==0.57143 & ...
                round(sensitivity_analysis_struct.sigma_perineurium,6)==round(8.70322e-4,6);
            color_data = sensitivity_analysis_struct.sigma_fill;
            color_data_label = '\sigma_{surround} (S/m)';
        case 2
            % Waveform vs. perineurial conductivity
            x_data = sensitivity_analysis_struct.common_time_vector_ms';
            x_data_label = 'time (ms)';
            subset_indices = round(sensitivity_analysis_struct.anisotropy_ratio_vals,2)==3.43 & ...
                round(sensitivity_analysis_struct.sigma_endoneurium_z,5)==0.57143 & ...
                round(sensitivity_analysis_struct.sigma_fill,6)==round(1.76,6);
            color_data = sensitivity_analysis_struct.sigma_perineurium;
            color_data_label = '\sigma_{perineurium} (S/m)';
        case 3
            % Waveform vs. longitudinal endoneural conductivity
            x_data = sensitivity_analysis_struct.common_time_vector_ms';
            x_data_label = 'time (ms)';
            subset_indices = round(sensitivity_analysis_struct.anisotropy_ratio_vals,2)==3.43 & ...
                round(sensitivity_analysis_struct.sigma_perineurium,6)==round(8.70322e-4,6) & ...
                round(sensitivity_analysis_struct.sigma_fill,6)==round(1.76,6);
            color_data = sensitivity_analysis_struct.sigma_endoneurium_z;
            color_data_label = '\sigma_z (S/m)';
        case 4
            % Waveform vs. anisotropy
            x_data = sensitivity_analysis_struct.common_time_vector_ms';
            x_data_label = 'time (ms)';
            subset_indices = round(sensitivity_analysis_struct.sigma_endoneurium_z,5)==0.57143 & ...
                round(sensitivity_analysis_struct.sigma_perineurium,6)==round(8.70322e-4,6) & ...
                round(sensitivity_analysis_struct.sigma_fill,6)==round(1.76,6);
            color_data = sensitivity_analysis_struct.anisotropy_ratio_vals;
            color_data_label = '\sigma_z / \sigma_r';
        otherwise
            error('Invalid case in switch statement');
    end
    
    y_data = ((sensitivity_analysis_struct.CAP_signals_all./max(abs(sensitivity_analysis_struct.CAP_signals_all)))');
    
    
    clear g
    g=gramm('x',x_data,...
        'y',y_data,...
        'color',color_data,...
        'subset',subset_indices);
    
    g.facet_grid([],round(360-sensitivity_analysis_struct.Theta_contact_CylUm300t_20230323_001,2));
    g.geom_line();
    
    g.set_names('column',['cuff opening (',char(176),')'],'x',x_data_label,'y',{'signal','(normalized)'},...
        'color',color_data_label);
    g.axe_property('XLim',x_bounds,'Visible','off');
    g.set_layout_options('legend',true,'legend_position','auto');
    g.set_text_options('interpreter','tex','base_size',12);
    %     g.set_color_options('map','parula')
    g.set_continuous_color('colormap','parula');
    figure('position',[680   558   612   420]);
    g.draw();
end


function [output_table, sensitivity_analysis_table] = load_sensitivity_analysis_results(flag_myelinated)

output_table = [];
sensitivity_analysis_table = [];

if (flag_myelinated==1)
    output_table_files = {
        '../results/mimic_ground_truth_myel_sweep_workspace.mat'
        };
elseif (flag_myelinated==0)
    output_table_files = {
        '../results/mimic_ground_truth_unmyel_sweep_workspace.mat'
        };
else
    error('invalid flag_myelinated value')
end


for i = 1:length(output_table_files)
    % Load the sensitivity analysis data and the sensitivity analysis table (with comma delimiter)
    if (flag_myelinated==1)
        workspace_i = load(output_table_files{i},'output_table');
        output_table_i = workspace_i.output_table;
    else
        workspace_i = load(output_table_files{i},'all_output_tables');
        output_table_i = vertcat(workspace_i.all_output_tables{:});
    end
    
    % Append the output table and the sensitivity analysis table
    output_table = [output_table;output_table_i];
end
