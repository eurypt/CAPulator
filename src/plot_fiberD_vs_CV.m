%{
    Plot Iint thresh vs. fiberD; fit equation
%}

% load('plot_fiberD_vs_SFAPs_output_from_20220301.mat')
myel_data = load('../bin/template_data_20221004_0_2022100409.mat');
unmyel_data = load('../bin/template_data_20221004_0_2022100410.mat');

all_fiber_types = [repmat({'myel'},length(myel_data.output_data_structure),1);repmat({'unmyel'},length(unmyel_data.output_data_structure),1)]';
all_fiber_diameters = [vertcat(myel_data.output_data_structure(:).fiber_diameter);vertcat(unmyel_data.output_data_structure(:).fiber_diameter)]';
CV_all = [vertcat(myel_data.output_data_structure(:).conduction_velocity_m_per_s);vertcat(unmyel_data.output_data_structure(:).conduction_velocity_m_per_s)]';

% myel
myel_indices=find(strcmp(all_fiber_types,'myel'));
x = all_fiber_diameters(myel_indices)';
y = CV_all(myel_indices)';

f_myel=fit(x,y,'poly1')
figure('color',[1 1 1]); plot(x,y,'k.'); %'LineWidth',2,'MarkerSize',20)
hold on;
x_eval = linspace(x(1),x(end),100);
plot(x_eval,feval(f_myel,x_eval),'k-');
% shadedErrorBar(x_eval,feval(f_myel,x_eval),0.05*feval(f_myel,x_eval));
xlabel('fiber diameter (\mum)');
ylabel('CV (m/s)');
title('myel')
legend({'data','fit'},'location','southeast')
set(gca,'FontSize',14)

% % unmyel
% unmyel_indices=find(strcmp(all_fiber_types,'unmyel'));
% x = all_fiber_diameters(unmyel_indices)';
% transformation_exponent = 2; % if set to 1, no transformation occurs
% if (transformation_exponent~=1)
%     warning('transformation_exponent~=1');
% end
% transformation = @(input) input.^(transformation_exponent);
% inv_transformation = @(input) input.^(1/transformation_exponent);
% y = CV_all(unmyel_indices)';
% y_transformed = transformation(y);
% 
% % x(1) = []; y_transformed(1) = []; warning('manually excluding data point index = 1 because its threshold is weird')
% % X_0 = 0.105287; % [um] manually specified smallest fiber diameter in the nerve (for verifying extrapolation of fit since initial sims at the lowets fiber diameter didn't yield valid Iint)
% 
% f_unmyel=fit(x,y_transformed,'poly1')
% figure('color',[1 1 1]); plot(x,inv_transformation(y_transformed),'k.'); %'LineWidth',2,'MarkerSize',20)
% hold on;
% x_eval = linspace(x(1),x(end),100);
% % shadedErrorBar(x_eval,inv_transformation(feval(f_unmyel,x_eval)),0.05*inv_transformation(feval(f_unmyel,x_eval)));
% plot(x_eval,inv_transformation(feval(f_unmyel,x_eval)),'k-');
% xlabel('fiber diameter (\mum)');
% ylabel('CV (m/s)');
% title('unmyel')
% legend({'data','fit'},'location','southeast')
% set(gca,'FontSize',14)
% if (transformation_exponent~=1)
%     warning('warning: transformation\_exponent~=1');
% end
% ylim_i = get(gca,'YLim');
% ylim([0 ylim_i(2)]); % set the bottom y axis value to zero if not already


% unmyel (transform x instead)
unmyel_indices=find(strcmp(all_fiber_types,'unmyel'));
x = all_fiber_diameters(unmyel_indices)';
transformation_exponent = 0.5; % if set to 1, no transformation occurs
if (transformation_exponent~=1)
    warning('transformation_exponent~=1');
end
transformation = @(input) input.^(transformation_exponent);
inv_transformation = @(input) input.^(1/transformation_exponent);
y = CV_all(unmyel_indices)';
x_transformed = transformation(x);

% x(1) = []; y_transformed(1) = []; warning('manually excluding data point index = 1 because its threshold is weird')
% X_0 = 0.105287; % [um] manually specified smallest fiber diameter in the nerve (for verifying extrapolation of fit since initial sims at the lowets fiber diameter didn't yield valid Iint)

f_unmyel_sqrtD=fit(x_transformed,y,'poly1')
figure('color',[1 1 1]); plot(inv_transformation(x_transformed),y,'k.'); %'LineWidth',2,'MarkerSize',20)
hold on;
x_eval = linspace(x(1),x(end),100);

% plot(x_eval,f_unmyel.p1*transformation(x_eval),'k-'); % plot wihout intercept
plot(x_eval,feval(f_unmyel_sqrtD,transformation(x_eval)),'k-');
xlabel('fiber diameter (\mum)');
ylabel('CV (m/s)');
title('unmyel')
legend({'data','fit'},'location','southeast')
set(gca,'FontSize',14)
ylim_i = get(gca,'YLim');
ylim([0 ylim_i(2)]); % set the bottom y axis value to zero if not already

% 
% %%
% 
% %{
%  calculate tstop needed based on CV for myel & unmyel
% %}
% 
% %%% Myel
% naive_tstop = 65;
% ms_offset = 10; % [ms] extra time to add on top of the length-based time; prevents super short sims for the fastest fibers
% percent_adjustment = 1.15; 
% ax_length = 45;
% adaptive_t_total_myel = round(sum((percent_adjustment*ax_length)./CV_all(myel_indices)) + ...
%     ms_offset*length(myel_indices)) % tstop adapats to a percentage of the length
% naive_t_total_myel = naive_tstop*length(CV_all(myel_indices)) % tstop = 120 [ms] for all fibers in x
% 
% %%% Unmyel
% naive_tstop = 120;
% adaptive_t_total_unmyel = round(sum((percent_adjustment*ax_length)./CV_all(unmyel_indices)) + ...
%     +ms_offset*length(myel_indices)) % tstop adapats to a percentage of the length
% naive_t_total_unmyel = naive_tstop*length(CV_all(unmyel_indices)) % tstop = 120 [ms] for all fibers in x
% 
% 
% %%
% 
% %{
% Calculate time between nodes (myel) or compartments (unmyel) based on CV
% and z coords
% %}
% 
% % myel
% myel_indices=find(strcmp(all_fiber_types,'myel'));
% n_compartments_per_repeat = 11; % assume that there are this many compartments for each node (i.e., total compartment counts of node, MYSA, FLUT, STIN)
% delta_z_um = cellfun(@(x) mode(diff(x(3,1:n_compartments_per_repeat:end))),coords_all(myel_indices))'; % [um] space between nodes (myel) or compartments (unmyel)
% x = all_fiber_diameters(myel_indices)';
% y = delta_z_um./CV_all(myel_indices)';
% 
% f_myel=fit(x,y,'poly9')
% figure('color',[1 1 1]); plot(x,y,'.'); %'LineWidth',2,'MarkerSize',20)
% hold on;
% x_eval = linspace(x(1),x(end),100);
% shadedErrorBar(x_eval,feval(f_myel,x_eval),0.05*feval(f_myel,x_eval));
% xlabel('fiber diameter (\mum)');
% ylabel({'time between nodes (myel)','or compartments (unmyel) (\mus)'});
% title('myel')
% legend({'data','fit'},'location','southeast')
% set(gca,'FontSize',14)
% 
% % unmyel
% unmyel_indices=find(strcmp(all_fiber_types,'unmyel'));
% n_compartments_per_repeat = 1; % assume that there are this many compartments for each node (i.e., total compartment counts of node, MYSA, FLUT, STIN)
% delta_z_um = cellfun(@(x) mode(diff(x(3,1:n_compartments_per_repeat:end))),coords_all(unmyel_indices))'; % [um] space between nodes (myel) or compartments (unmyel)
% x = all_fiber_diameters(unmyel_indices)';
% transformation_exponent = 1; % if set to 1, no transformation occurs
% if (transformation_exponent~=1)
%     warning('transformation_exponent~=1');
% end
% transformation = @(input) input.^(transformation_exponent);
% inv_transformation = @(input) input.^(1/transformation_exponent);
% y = delta_z_um./CV_all(unmyel_indices)';
% y_transformed = transformation(y);
% 
% % x(1) = []; y_transformed(1) = []; warning('manually excluding data point index = 1 because its threshold is weird')
% % X_0 = 0.105287; % [um] manually specified smallest fiber diameter in the nerve (for verifying extrapolation of fit since initial sims at the lowets fiber diameter didn't yield valid Iint)
% 
% f_unmyel=fit(x,y_transformed,'poly9')
% figure('color',[1 1 1]); plot(x,inv_transformation(y_transformed),'.'); %'LineWidth',2,'MarkerSize',20)
% hold on;
% x_eval = linspace(x(1),x(end),100);
% shadedErrorBar(x_eval,inv_transformation(feval(f_unmyel,x_eval)),0.05*inv_transformation(feval(f_unmyel,x_eval)));
% xlabel('fiber diameter (\mum)');
% ylabel({'time between nodes (myel)','or compartments (unmyel) (\mus)'});
% title('unmyel')
% legend({'data','fit'},'location','southeast')
% set(gca,'FontSize',14)
% if (transformation_exponent~=1)
%     title('warning: transformation\_exponent~=1');
% end
% 
