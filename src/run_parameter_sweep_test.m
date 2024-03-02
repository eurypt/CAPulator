
function run_parameter_sweep_test(params_json_filename,parameter_sweep_values,...
    parameter_sweep_variable_name,parameter_sweep_xlabel,title_labels,reference_CNAP_label,x_bounds)
% Load the base JSON to be tweaked durign this sensitivity analysis
params = loadjson(params_json_filename);

V_pk2pk = [];
for i = 1:length(parameter_sweep_values)
    params.(parameter_sweep_variable_name) = parameter_sweep_values(i);
    output_table = process_params(params);
    
    % calculate the peak-to-peak voltage
    % Load reference signal
    % ref = load(reference_signals{i},'time_vector','CNAP_signal_uV');
    match_ind = find(strcmp(output_table.label,reference_CNAP_label));
    ref.time_vector = output_table.common_time_vector_ms{match_ind};
    ref.CNAP_signal_uV = output_table.CAP_signal_uV{match_ind};
    reference_signal_time_vector = ref.('time_vector');
    % Ensure reference_signal_time_vector is a column vector
    reference_signal_time_vector = reference_signal_time_vector(:);
    reference_signal = ref.('CNAP_signal_uV');
    V_pk2pk = [V_pk2pk,nan(size(output_table,1),1)];
    %assert(~isempty(reference_ind) && length(reference_ind)==1, 'invalid value of reference_ind; it should be a non-empty scalar');
    for j = 1:size(output_table,1)
        
        % calculate the rmse between the reference and the present CNAP
        % use interp1 mostly to zero-pad or truncate appropriately so the
        % test signal is the same length as the reference signal
        test_signal = interp1(output_table.common_time_vector_ms{j},output_table.CAP_signal_uV{j},...
            reference_signal_time_vector,'linear',0);
        V_pk2pk(j,i) = 100*max(abs(test_signal-reference_signal))/(max(reference_signal)-...
            min(reference_signal));
        %test_signal = output_table.CAP_signal_uV{j};
        %assert(mode(diff(output_table.common_time_vector_ms{j}))==mode(diff(output_table.common_time_vector_ms{reference_ind})));
        %rmse_array(j) = sqrt(mean((reference_signal(test_time_indices)-test_signal(test_time_indices)).^2));
        %{
            V_pk2pk(j,i) = max(test_signal)-min(test_signal);
        %}
    end
    
    % Plot all the output_table signal overlaid
    fig1 = custom_figure;
    for j = 1:size(output_table,1)
        custom_plot(output_table.common_time_vector_ms{j},output_table.CAP_signal_uV{j},'time (ms)','CNAP signal (\muV)');
        hold on;
    end
    % % Plot the ground truth signal on top
    % plot(reference_signal_time_vector,reference_signal,'k-')
    % legend([output_table.label;reference_signals],'location','best')
    legend(output_table.label,'location','best')
    title(title_labels{i})
    % Make the ylim symmetric
    ylim_i = get(gca,'YLim'); 
    ylim(max(abs(ylim_i))*[-1 1])
    % Set x bounds 
    if ~isempty(x_bounds{i})
        xlim(x_bounds{i});
    end

    % Save figure
    [~,name,~] = fileparts(params_json_filename);
    fig_i_output_name = ['../results/',name,'_',strrep(strrep(title_labels{i},'.','p'),' ','_')];
    warning('figure name is not necessarily unique: %s',fig_i_output_name)
    savefig(fig_i_output_name)
    
    % Close then open the figure so that if you edit manually you can
    % easily save the edits by clicking save rather than having to match up
    % the figure with figure name
    close(fig1);
    openfig(fig_i_output_name);
end


% plot results
fig2 = custom_figure('position',[625,515,979,416]);
custom_plot(1:size(V_pk2pk,1),V_pk2pk,'.-','MarkerSize',20,'',...
    {'Max % Discrepancy',['Relative to CNAP_{',reference_CNAP_label,'}']})

tick_vals =  4.^[-1:6]; 
tick_labels = arrayfun(@(x) num2str(x,'%0.3f'),tick_vals,'UniformOutput',false); 
set(gca,'YScale','log','YTick',tick_vals,'YTickLabel',tick_labels);
set(gca,'XTick',1:size(output_table,1),'XTickLabel',output_table.label(1:end))
legend(title_labels(1:size(V_pk2pk,2)),'location','best')
xlabel(parameter_sweep_xlabel)
set(gcf,'Position',[108   236   494   336])
% colororder(parula(size(V_pk2pk,2)));
hold on;
% Save summary figure
summary_fig_output_name = ['../results/',name,'_sweep_summary.fig'];
warning('figure name is not necessarily unique: %s',summary_fig_output_name)
savefig(summary_fig_output_name);

% Close then open the figure so that if you edit manually you can
% easily save the edits by clicking save rather than having to match up
% the figure with figure name
close(fig2)
openfig(summary_fig_output_name)
