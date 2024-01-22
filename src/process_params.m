%{

Calculate the parameter sweep for the CNAP models based on the present of param_sweep_group in the params struct
If param_sweep_group is present, then the function will call itself recursively, creating different versions of params
If param_sweep_group is not present, then the function will call calculate_CAP_via_temporal_templates to calculate the CAP signal

Inputs:
params: struct containing the parameters for the CAP calculation

Outputs:
output_table: table containing the results of the CAP calculation

%}
function output_table = process_params(params)

% If param_sweep_group is a field in params, create different versions of
% params (named params_i) that contain the different sets of field values
% present in the param_sweep_group; then use each of these versions of
% params to call the present function recursively
output_table = [];
if (isfield(params,'param_sweep_group'))
    for i = 1:length(params.param_sweep_group)
        params_i = rmfield(params,'param_sweep_group');
        % determine whether the param_sweep_group is a cell array or a
        % struct array; the former happens if the different elements of
        % params_sweep_group have different sets of fields (e.g., if one
        % element has fields x,y,z, while the another element only has
        % field z.
        if (iscell(params.param_sweep_group))
            param_sweep_group_i = params.param_sweep_group{i};
        elseif (isstruct(params.param_sweep_group))
            param_sweep_group_i = params.param_sweep_group(i);
        else
            error('unexpected format for params.param_sweep_group; expected it to be either a cell array or a struct array...')
        end
        f = fieldnames(param_sweep_group_i);
        for j = 1:length(f)
            assert(~isfield(params_i,f{j}),'Error: The following ''param_sweep_group'' field is also defined in a parent scope: %s',f{j});
            params_i.(f{j}) = param_sweep_group_i.(f{j});
        end
        if (isempty(output_table))
            output_table = process_params(params_i);
        else
            output_table = [output_table; process_params(params_i)];
        end
    end
else
    % [CAP_signal_uV, common_time_vector_ms, SFAPs_uV] = ...
    %     calculate_CAP_via_temporal_templates(params);
    [CAP_signal_uV, common_time_vector_ms, SFAPs_uV] = ...
        CAPulator.calculate_CAP_via_temporal_templates(params);
    
    output_struct = [];

    % if saving SFAPs was set to 0, delete SFAPs to save disk space
    if (isfield(params,'save_SFAPs') && params.save_SFAPs == 1)
        output_struct.SFAPs_uV = SFAPs_uV;
    end
    % save output to table
    % note: the label field is used to identify the parameters used to generate the signal
    if (isfield(params,'label'))
        output_struct.label = params.label;
    end

    output_struct.params = params;
    output_struct.CAP_signal_uV = CAP_signal_uV;
    output_struct.common_time_vector_ms = common_time_vector_ms;

    % Convert the output struct to a table
    output_table = struct2table(output_struct,'AsArray',true);
    
end
