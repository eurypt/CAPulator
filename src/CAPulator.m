%{
    This is a class that has only static methods. It is used to calculate the
    CNAP signal due to the specified params. The following methods are contained
    in this class:

    calculate_CAP_via_temporal_templates - calculates the CNAP signal due to the
        specified params
    get_template_data - gets the template data from the specified file
    get_straight_fiber_xyz_coords - gets the xyz coordinates of the specified
        fibers
    get_electric_potentials - gets the electric potentials at the specified
        locations due to the specified extracellular recording model
    
%}
classdef CAPulator
    methods(Static)
        %{
            This function is the main function of this class. It calculates the
            CNAP signal due to the specified params. The following are the inputs

            template_data_source_filename - the filename of the file containing the
                template data
            fiber_type - the type of fiber (myel or unmyel)
            xy_coords_um - the xy coordinates of the fibers
            z_min_um - the minimum z coordinate of the fibers
            z_max_um - the maximum z coordinate of the fibers
            extracellular_recording_model_filename - the filename of the file
                containing the extracellular recording model
            fiber_diameters_in_nerve_um - the fiber diameters in the nerve in um
            stim_location_change_mm - the change in stimulation location in mm
            xy_coords_filename - the filename of the file containing the xy
                coordinates of the fibers (optional, but either xy_coords_filename
                or xy_coords_um must be specified)
            xy_coords_um - the xy coordinates of the fibers (optional but either
                xy_coords_filename or xy_coords_um must be specified)

            The following are the optional inputs

            flag_dipole - if 1, use the dipole approach; if 0, use the monopole
                approach
            interp_type - the type of interpolation to use when interpolating the
                CVs of the fibers
            flag_jitter_z - if 1, jitter the z coordinates of the fibers
            diameter_precision - the precision to use when rounding the fiber
                diameters

                

        %}
        function [CNAP_signal_uV, final_time_vector_ms, SFAPs_uV] = ...
                calculate_CAP_via_temporal_templates(params)

            % Extract all the required parameter field values
            template_data_source_filename = params.template_data_source_filename;

            % fiber_type=params.fiber_type;
            extracellular_recording_model_filename = params.extracellular_recording_model_filename;

            LENGTH_AXONS_um = params.z_max_um-params.z_min_um;
            DESIRED_CENTER_um = mean([params.z_max_um,params.z_min_um]);

            fiber_diameters_in_nerve_um = params.fiber_diameters_in_nerve_um;

            % Load xy_coords_um from the file if it was specified; otherwise, use the
            % specified xy_coords_um
            if (isfield(params,'xy_coords_filename'))
                xy_coords_filename = params.xy_coords_filename;
                load(xy_coords_filename,'xy_coords_um');
            else
                assert(isfield(params,'xy_coords_um'),'xy_coords_um or xy_coords_filename must be specified');
                xy_coords_um=params.xy_coords_um;
            end

            stim_location_change_mm = params.stim_location_change_mm;

            % Done extracting all the required parameter field values

            % Extract all the optional parameter field values
            if (isfield(params,'flag_dipole'))
                flag_dipole=params.flag_dipole;
            else
                % default to dipole approach
                flag_dipole = 1;
            end

            if (isfield(params,'template_diameters_to_use_um'))
                template_diameters_to_use_um=params.template_diameters_to_use_um;
            else
                % default to using all template diameters
                template_diameters_to_use_um = [];
            end

            if (isfield(params,'interp_type'))
                interp_type=params.interp_type;
            else
                % default to using all template diameters
                interp_type = 'linear';
            end
            if (isfield(params,'flag_jitter_z'))
                flag_jitter_z = params.flag_jitter_z;
            else
                flag_jitter_z = 0;
            end
            if (isfield(params,'dset_tag'))
                dset_tag = params.dset_tag;
            else
                dset_tag = 'dset1';
            end
            if (isfield(params,'outersolnum'))
                outersolnum = params.outersolnum;
            else
                outersolnum = 1;
            end
            if (isfield(params,'diameter_precision'))
                diameter_precision = params.diameter_precision;
            else
                diameter_precision = 3;
            end

            % Done extracting all the optional parameter field values


            % ensure fiber_diameters_in_nerve_um is a row vector
            if (size(fiber_diameters_in_nerve_um,1)>1)
                fiber_diameters_in_nerve_um = transpose(fiber_diameters_in_nerve_um);
            end

            % Round fiber diameters to specified precision
            fiber_diameters_in_nerve_um = round(fiber_diameters_in_nerve_um,diameter_precision);

            % transpose xy_coords_um so that it is a 2-by-M matrix in which M is the
            % number of fibers *or* in which M=1 (i.e., one xy coordinate for all
            % axons)
            xy_coords_um = transpose(xy_coords_um);

            %{
                If the following fields are not already defined in the params
                structure, then get them by running the get_template_data method
                of this class on the specified template_data_source_filename.
                This method returns the following outputs:
                    'all_conduction_velocity_m_per_s','all_fiber_diameters',...
                    'common_time_vector_ms','template_data_matrix','reference_peak_times_all_templates',...
                    'reference_compartment_locations_um_all_templates','fiber_type'
                
                If those fields are already defined in the params structure,
                then just use them (define the variables outside of the
                structure)
            %}
            if (~isfield(params,'all_conduction_velocity_m_per_s') || ...
                    ~isfield(params,'all_fiber_diameters') || ...
                    ~isfield(params,'common_time_vector_ms') || ...
                    ~isfield(params,'template_data_matrix') || ...
                    ~isfield(params,'reference_peak_times_all_templates') || ...
                    ~isfield(params,'reference_compartment_locations_um_all_templates') || ...
                    ~isfield(params,'fiber_type'))
                [all_conduction_velocity_m_per_s,all_fiber_diameters,...
                    common_time_vector_ms,template_data_matrix,reference_peak_times_all_templates,...
                    reference_compartment_locations_um_all_templates,fiber_type] = ...
                    CAPulator.get_template_data(template_data_source_filename,flag_dipole);
            else
                all_conduction_velocity_m_per_s = params.all_conduction_velocity_m_per_s;
                all_fiber_diameters = params.all_fiber_diameters;
                common_time_vector_ms = params.common_time_vector_ms;
                template_data_matrix = params.template_data_matrix;
                reference_peak_times_all_templates = params.reference_peak_times_all_templates;
                reference_compartment_locations_um_all_templates = params.reference_compartment_locations_um_all_templates;
                fiber_type = params.fiber_type;
            end

            %

            % Downsample the templates to just above 2x the nyquist frequency as determined
            % by precondition_templates.m. Based on that, the nyquist frequency is 47.1 kHz
            % for myelinated fibers and 7.8kHz for unmyelinatd fibers. So, cut out all
            % frequencies above 50 kHz and 8 kHz (via spectral subtraction) then downsample
            % to the signals to 16 kHz for unmyelinated fibers and 100 kHz for unmyelinated
            % fibers.
            if (strcmp(fiber_type,'unmyelinated'))
                [fft_i,freq] = CAPulator.fft_wrapper(template_data_matrix,[],[],1/(1e-3*mode(diff(common_time_vector_ms))));
                cutoff_freq_Hz = 8e3;
                fft_i2 = fft_i(abs(freq)<=cutoff_freq_Hz,:);
                amplitude_adjustment = size(fft_i2,1)/size(template_data_matrix,1);
                template_data_matrix = ifft(fft_i2)*amplitude_adjustment;
                new_sampling_rate_ms = (1e3*(1/cutoff_freq_Hz)/2); % [ms]
                common_time_vector_ms = common_time_vector_ms(1) + ((1:size(template_data_matrix,1))-1)*new_sampling_rate_ms;
            end

            % Round all_fiber_diameters to the specified precision
            all_fiber_diameters = round(all_fiber_diameters,diameter_precision);

            % manually declare the variables as below because parfor may not be able to
            % access variables loaded directly from a file
            % Use only the specified subset of fiber diameters as templates
            if (isempty(template_diameters_to_use_um))
                template_diameters_to_use_um = unique(all_fiber_diameters);
            end
            templates_to_keep = find(ismember(all_fiber_diameters,...
                round(template_diameters_to_use_um,diameter_precision)));
            template_data_matrix = template_data_matrix(:,templates_to_keep,:);

            all_fiber_diameters = all_fiber_diameters(templates_to_keep);
            all_conduction_velocity_m_per_s = all_conduction_velocity_m_per_s(templates_to_keep);
            reference_peak_times_all_templates = reference_peak_times_all_templates(templates_to_keep);
            reference_compartment_locations_um_all_templates = reference_compartment_locations_um_all_templates(templates_to_keep);
            % common_time_vector_ms = common_time_vector_ms;

            % assert(length(expected_initiation_delay_ms)==1 || length(expected_initiation_delay_ms)==size(template_data_matrix,2),...
            %     'expected_initiation_delay_ms should contain either one value per template*or* a single scalar value');
            assert(~isempty(all_fiber_diameters),'all_fiber_diameters should not be empty; verify that template_data_source_filename is set appropriately');


            % Remove any fiber diameters that exceed the bounds of the template
            % transmembrane currents
            diameters_to_remove = find(fiber_diameters_in_nerve_um<min(all_fiber_diameters) | ...
                fiber_diameters_in_nerve_um>max(all_fiber_diameters));
            if (~isempty(diameters_to_remove))
                warning('removing %d fibers because they are outside the bounds of the fiber diameters of available Im',...
                    length(diameters_to_remove));
                fiber_diameters_in_nerve_um(diameters_to_remove) = [];

                % Remove xy coordinates of the removed fiber if there were multiple xy
                % coordinates provided
                if (size(xy_coords_um,2)>1)
                    xy_coords_um(diameters_to_remove) = [];
                end
            end

            query_fiber_diameter_interpolated_all = interp1(all_fiber_diameters,all_fiber_diameters,...
                fiber_diameters_in_nerve_um,interp_type);
            [unique_fiber_diameters,n_of_each_unique_fiber_diameter,IA,IC] = ...
                CAPulator.count_unique_values(query_fiber_diameter_interpolated_all,diameter_precision);


            % Interpolate CV
            % If the fiber type is myelinated, interpolate CV linearly across fiber
            % diameters; if the fiber type is unmyelinated, interpolate CV *squared* linearly
            % across fiber diameters. Since the square of the CV for unmyelinated fibers is a
            % linear function of fiber diameter, this approach produces better CV estimates
            if (strcmp(fiber_type,'myelinated'))
                CV_interpolated_all_fibers = interp1(all_fiber_diameters,...
                    all_conduction_velocity_m_per_s,...
                    query_fiber_diameter_interpolated_all,interp_type);
            elseif (strcmp(fiber_type,'unmyelinated'))
                CV_interpolated_all_fibers = sqrt(interp1(all_fiber_diameters,...
                    all_conduction_velocity_m_per_s.^2,...
                    query_fiber_diameter_interpolated_all,interp_type));
            else
                error('fiber_type must be either myelinated or unmyelinated');
            end
            CV_interpolated_at_unique_diameters = CV_interpolated_all_fibers(IA);

            % ===== Start this major output: XYZ coordinates of all fibers =====
            % Get xyz coordinates of all compartments; since the temporal template
            % method depends on repeating a template for each node, and since there can
            % be multiple compartments per node (namely in myelinated fibers), truncate
            % the list of coordinates so that its length is an integer multiple of
            % n_compartments_per_node

            % If coords_all_compartments_all_axons, n_compartments_per_node, and
            % distance_between_nodes_all_fibers_um are already fields in the params, then use
            % those values; otherwise, calculate them
            if (isfield(params,'coords_all_compartments_all_axons') && ...
                    isfield(params,'n_compartments_per_node') && ...
                    isfield(params,'distance_between_nodes_all_fibers_um'))
                coords_all_compartments_all_axons = params.coords_all_compartments_all_axons;
                n_compartments_per_node = params.n_compartments_per_node;
                distance_between_nodes_all_fibers_um = params.distance_between_nodes_all_fibers_um;
            else
                geometry_determination_method = 2; % geometry_determination_method = 0 for preset fiber diameters; geometry_determination_method = 1 for MRG-based geometry interpolation; geometry_determination_method = 2 for GeometryBuilder fits from SPARC Y2Q1
                [coords_all_compartments_all_axons,n_compartments_per_node,distance_between_nodes_all_fibers_um] = ...
                    CAPulator.get_straight_fiber_xyz_coords(fiber_type,geometry_determination_method,...
                    query_fiber_diameter_interpolated_all,xy_coords_um,LENGTH_AXONS_um,DESIRED_CENTER_um,[],flag_jitter_z);
            end
            % ===== End this Major Output =====

            distance_um_between_nodes_unique_fiber_diameters = distance_between_nodes_all_fibers_um(IA);

            time_between_nodes_ms = 1e-3*distance_between_nodes_all_fibers_um'./CV_interpolated_all_fibers;

            % specify the volume conductor information to extract electric
            % potentials
            if (strcmp(extracellular_recording_model_filename(end-3:end),'.mat'))
                % If a .mat file was specified, it is assumed that all
                % fibers are at a single xy coordinate (e.g., all fibers
                % placed at the origin); then the potentials will be
                % interpolated from the z_coord and V in the .mat file
                extracellular_recording_model_data = load(extracellular_recording_model_filename,'z_coords','V');
                extracellular_recording_model_structure.z_coords = extracellular_recording_model_data.z_coords;
                extracellular_recording_model_structure.V = extracellular_recording_model_data.V;
            elseif (strcmp(extracellular_recording_model_filename(end-3:end),'.mph'))
                extracellular_recording_model_structure.filename = extracellular_recording_model_filename;
                extracellular_recording_model_structure.outersolnum = outersolnum;
                extracellular_recording_model_structure.dset_tag = dset_tag;
            else
                error('extracellular_recording_model_filename should be either a .mat file (containing z_coords & V in a straight line) *or* a COMSOL (.mph) file');
            end
            % get the electric potentials using the volume conductor
            % information specified above
            [recording_sensivitity_all_fibers,n_compartments_per_fiber] = CAPulator.get_electric_potentials(...
                extracellular_recording_model_structure,coords_all_compartments_all_axons);

            % Reshape & permute the forward model to store the other repeating compartments
            % (e.g., MYSA, FLUT, STIN) in the third dimension.
            % Also, if using the dipole approach, iterate through every element in
            % recording_sensivitity_all_fibers to calculate the diff
            for i = 1:length(recording_sensivitity_all_fibers)
                if (flag_dipole)
                    recording_sensivitity_all_fibers{i} = [diff(recording_sensivitity_all_fibers{i});0];
                end
                recording_sensivitity_all_fibers{i} = reshape(recording_sensivitity_all_fibers{i},n_compartments_per_node,[]);
                recording_sensivitity_all_fibers{i} = permute(recording_sensivitity_all_fibers{i},[2,1]);
            end


            % Pre-identify all the closest template fiber diameters (upper and lower)
            % for each query fiber diameter and weights; this approach elegantly
            % handles when fiber diameters match a template since one of the weights is
            % set to zero and the other is set to one
            fiber_floating_point_indices = interp1(all_fiber_diameters,1:length(all_fiber_diameters),...
                query_fiber_diameter_interpolated_all',interp_type);
            left_integer_indices_and_weights = [floor(fiber_floating_point_indices),1-rem(...
                fiber_floating_point_indices,1)];
            right_integer_indices_and_weights = [ceil(fiber_floating_point_indices),rem(...
                fiber_floating_point_indices,1)];
            left_integer_indices_and_weights_unique_diams = left_integer_indices_and_weights(IA,:);
            right_integer_indices_and_weights_unique_diams = right_integer_indices_and_weights(IA,:);

            %%% Iterate through each fiber diameter to get SFAP signals for each one
            dt_ms_original = mode(diff(common_time_vector_ms));



            % Calculate delays in terms of
            delays_ms_forward_model = time_between_nodes_ms.*(n_compartments_per_fiber/n_compartments_per_node);

            % length to expand templates to so that they can handle shifts due to
            % forward model
            new_length = ceil(size(template_data_matrix,1)+max(delays_ms_forward_model)/dt_ms_original);

            sampling_rate_Hz = 1e3/dt_ms_original; % [Hz]
            [template_data_matrix_fft,template_freq_Hz] = CAPulator.fft_wrapper(template_data_matrix,...
                new_length,[],sampling_rate_Hz);

            SFAP_all_fibers_fft = 1i*ones(size(template_data_matrix_fft,1),length(recording_sensivitity_all_fibers));

            all_SFAP_delays_ms = zeros(length(recording_sensivitity_all_fibers),1);
            for unique_fiber_diameter_ind = length(unique_fiber_diameters):-1:1

                % Define the left weights and right weights and indinces
                left_ind = left_integer_indices_and_weights_unique_diams(unique_fiber_diameter_ind,1);
                left_weight = left_integer_indices_and_weights_unique_diams(unique_fiber_diameter_ind,2);
                right_ind = right_integer_indices_and_weights_unique_diams(unique_fiber_diameter_ind,1);
                right_weight = right_integer_indices_and_weights_unique_diams(unique_fiber_diameter_ind,2);

                % Combine forward models for fibers of the same fiber diameter This assumes
                % that fibers of the same fiber diameter have the same number of
                % compartments and that activation happens at the same location in those
                % fibers
                indices_to_combine = find(round(unique_fiber_diameters(unique_fiber_diameter_ind),diameter_precision)==...
                    round(query_fiber_diameter_interpolated_all,diameter_precision));
                % pull out the z coordinates of the present unique fiber diameter; this
                % assumes that all the axons of a given fiber diameter are straight and
                % of equal number of compartments, but necessarily of the same xy
                % coordinates
                nodes_z_coords_um = coords_all_compartments_all_axons{indices_to_combine(1)}(3,1:n_compartments_per_node:end);
                assert(all(cellfun(@(x) size(x,2)==size(coords_all_compartments_all_axons{indices_to_combine(1)},2),...
                    coords_all_compartments_all_axons(indices_to_combine))),...
                    ['currently the way that the index_of_node_where_first_AP_occurred is ',...
                    'set assumes all fibers of a given unique fiber diameter have the same ',...
                    'number of compartments and are straight fibers']);

                % prep a variable for the forward model; also count number of
                % nodes; confirm that the length of the forward model is an integer
                % multiple of n_compartments_per_node since the temporal template gets
                % repeated for each node
                %     forward_model_i = recording_sensivitity_all_fibers{fiber_ind};
                %     forward_model_i = sum(cell2mat(transpose(recording_sensivitity_all_fibers(indices_to_combine))),2);
                % Keep individual fibers separate such that the 2D matrix has columns
                % of one fiber followed by columns of the next fiber, etc.
                forward_model_i = cell2mat(transpose(recording_sensivitity_all_fibers(indices_to_combine)));

                time_between_nodes_ms_i = unique(time_between_nodes_ms(indices_to_combine));
                assert(length(time_between_nodes_ms_i)==1,'The time between nodes should be the same for all fibers of the same diameter...');
                desired_dt_resampled = time_between_nodes_ms_i / ceil(time_between_nodes_ms_i/dt_ms_original);
                dt_ms = desired_dt_resampled;

                % What is the peak time at the first modeled node? If it is positive,
                % (a) make a delay vector from 0 to the # of nodes, (b) multiply it by the
                % conduction delay samples between nodes, and then (c) shift that whole
                % delay vector by round(peak time at first modeled node / time between
                % nodes).
                % If it is negative, calculate index of the first node in which the
                % delay is positive, then take the absolute value of that node index
                % subtracted from the 0 to # of nodes delay vector and proceed to (b)
                % from above
                % fit the inverse of peak times since the peak time is inversly
                % proportional to CV, and CV is linearly proportional to fiber
                % diameter. This gives a noticeably better estimate of peak time at
                % small fiber diameters when using low template counts while also
                % having no clear negative impact at the larger fiber diameters
                %     interp_peak_time_ms = 1/(left_weight*1/reference_peak_times_all_templates(left_ind) + ...
                %         right_weight*1/reference_peak_times_all_templates(right_ind)); % [ms]
                if (strcmp(fiber_type,'myelinated'))
                    interp_peak_time_ms = 1/interp1(all_fiber_diameters,1./reference_peak_times_all_templates,...
                        unique_fiber_diameters(unique_fiber_diameter_ind));

                elseif (strcmp(fiber_type,'unmyelinated'))
                    % Assuming the starting location and the ref location are the same for
                    % all axons (i.e., delta_x is constant), and knowing that CV *squared*
                    % for unmyelianted fibers is a linear function of fiber diameter, the
                    % best estimate of interp_peak_time_ms for unmyelinated fibers is
                    % obtained by linearly interpolating between the square of the
                    % reciprocal of the known points.
                    interp_peak_time_ms = sqrt(1/interp1(all_fiber_diameters,(1./reference_peak_times_all_templates).^2,...
                        unique_fiber_diameters(unique_fiber_diameter_ind)));

                end

                % Adjust the interp peak time by the time it takes to traverse the
                % change in stimulation distance; a positive shift in location moves
                % the stimulation closer to the recording, thus reducing the interp
                % peak time
                interp_peak_time_ms = interp_peak_time_ms - stim_location_change_mm / CV_interpolated_at_unique_diameters(unique_fiber_diameter_ind);
                approximate_ref_node_loc_um = (left_weight*reference_compartment_locations_um_all_templates(left_ind) + ...
                    right_weight*reference_compartment_locations_um_all_templates(right_ind));
                % increment by one since a difference of 0 indicates that node #1 is
                % the ref node
                interp_ref_node_ind = 1 + round((approximate_ref_node_loc_um - nodes_z_coords_um(1))/...
                    distance_um_between_nodes_unique_fiber_diameters(unique_fiber_diameter_ind));

                % substract 1 since if interp_ref_node_ind is equal to ref node ind
                % then no need to adjust interp_peak_ms
                peak_time_at_first_modeled_node_ms = interp_peak_time_ms - time_between_nodes_ms_i*(interp_ref_node_ind-1);

                % start the delay vector (in terms of nodes) from 0 to (# nodes modeled - 1)
                %     number_of_nodes_to_model = length(forward_model_i)/n_compartments_per_node;
                number_of_nodes_to_model = size(forward_model_i,1);

                % set a bit arbitrarily based on the assumption that the signal was
                % resampled to have a time between nodes that is an integer multiple of
                % sampling rate
                desired_dt_resampled = time_between_nodes_ms_i / ceil(time_between_nodes_ms_i/dt_ms_original);
                ROUND_PRECISION_FOR_CONDUCTION_SAMPLES_PER_NODE = 3;
                samples_between_nodes = round(time_between_nodes_ms_i/desired_dt_resampled,...
                    ROUND_PRECISION_FOR_CONDUCTION_SAMPLES_PER_NODE);

                % Store all the delays for each SFAP (in ms)
                % all_SFAP; adjust for the template initial time (i.e.,
                % common_time_vector_ms(1))
                all_SFAP_delays_ms(indices_to_combine) = peak_time_at_first_modeled_node_ms;

                % if the conduction delay between compartments/nodes is an
                % approximately an integer multiple of dt, then use a sparse shift
                % matrix of ones and zeros since it is more efficient that the
                % genalized sinc approach; otherwise, use the sinc approach to
                % calculate shift matrix
                % add 1 so that a delay of 0 converts into an index of 1, etc.
                i_indices = 1 + samples_between_nodes*(0:(number_of_nodes_to_model-1));
                j_indices = 1:number_of_nodes_to_model;
                s = ones(size(i_indices));
                all_shifts = sparse(i_indices,j_indices,s);


                % reshape the forward model in preparation for weighting; as with the
                % number of nodes calculations above, exclude the very last compartment;
                % reshape so that the weighting can occur for myelinated fibers that have
                % multiple repeating compartments per node
                %     reshaped_forward_model = transpose(reshape(forward_model_i,n_compartments_per_node,[]));
                % weight the shifts, then use the result as a filter for the template; do
                % this separate for each compartment of the repeating unit
                weighted_all_shifts = all_shifts*forward_model_i;

                % reshape & permute forward model directly
                % temp = reshape(forward_model_i,size(forward_model_i,1),n_compartments_per_node,[]);
                temp = reshape(weighted_all_shifts,size(weighted_all_shifts,1),n_compartments_per_node,[]);
                temp2 = permute(temp,[1 3 2]);

                % set the nfft so that the frequency resolution equals that of the
                % fft of temporal templates which has a frequency resolution equal
                % to (1/dt_ms_original)/length(final_expected_time_vector)
                % when dt_ms=time_between_nodes_ms, the original forward model can
                % be used for the fft
                dt_ms = desired_dt_resampled; %time_between_nodes_ms_i; %ceil(time_between_nodes_ms_i*cutoff_freq_Hz);
                nfft = round((1/dt_ms)*size(template_data_matrix_fft,1)/(1/dt_ms_original));

                % take FFT of forward model (without the upsamples)
                [filter_fft,filter_freq_Hz] = CAPulator.fft_wrapper(temp2,nfft,[],1e3/dt_ms);

                % if the filter length is shorter than the templates, calculate padding
                % needed to add and add it in
                padding_length = size(template_data_matrix_fft,1)-size(filter_fft,1);
                % If padding needed is positive, pad the fft; otherwise, truncate it
                if (padding_length>0)
                    fft_midpoint = ceil((nfft+1)/2);
                    filter_fft_padded = [filter_fft(1:fft_midpoint,:,:); ...
                        zeros(padding_length,size(filter_fft,2),size(filter_fft,3)); ...
                        filter_fft((fft_midpoint+1):end,:,:)];
                    % if the pre-padded fft is even, adjust the values accordingly
                    if  mod(size(filter_fft,1),2)==0
                        filter_fft_padded(fft_midpoint,:,:) = filter_fft(fft_midpoint,:,:)/2;
                        filter_fft_padded(fft_midpoint+padding_length,:,:) = filter_fft(fft_midpoint,:,:)/2;
                    end
                elseif (padding_length<0)
                    N = size(template_data_matrix_fft,1);
                    filter_fft_padded = filter_fft(([1:ceil(N/2),(end-(floor(N/2)-1)):end]),:,:);
                    % if the pre-padded fft is even, adjust the values accordingly
                    if  mod(size(filter_fft_padded,1),2)==0
                        filter_fft_padded(ceil(N/2)+1,:,:) = filter_fft(ceil(N/2)+1,:,:)+filter_fft(end-(floor(N/2)-1),:,:);
                    end
                    % assert(isreal(ifft(filter_fft_padded)));
                else
                    filter_fft_padded = filter_fft;
                end

                temporal_template_i_fft = left_weight*template_data_matrix_fft(:,left_ind,:) + ...
                    right_weight*template_data_matrix_fft(:,right_ind,:);

                SFAP_all_fibers_fft(:,indices_to_combine) = ...
                    sum(filter_fft_padded.*temporal_template_i_fft,3);
            end


            SFAPs_uV = 1e6*ifft(SFAP_all_fibers_fft); % [uV]

            % Add the template initial time offset to the delay correction
            all_SFAP_delays_ms = all_SFAP_delays_ms + common_time_vector_ms(1);

            % Adjust the SFAP delays by a specified amount before converting back into
            % time domain
            % Specify shifts
            % assert(all(all_SFAP_delays_ms>0),'The code expects only positive delays');
            all_SFAP_delays_samples = round(all_SFAP_delays_ms/dt_ms_original);
            % Zero pad the end of the SFAPs by the maximum shift to accomodate the
            % shifts
            min_delay_samples = min(all_SFAP_delays_samples);
            if (min_delay_samples>0)
                min_delay_samples = 0;
            end
            max_delay_samples = max(all_SFAP_delays_samples);
            SFAPs_uV = [
                zeros(-min_delay_samples,size(SFAPs_uV,2))
                SFAPs_uV
                zeros(max_delay_samples,size(SFAPs_uV,2))
                ];
            % Iterate through all the fibers and shift each
            n_fibers = size(SFAPs_uV,2);
            for i = 1:n_fibers
                if (all_SFAP_delays_samples(i)<0)
                    SFAPs_uV(:,i) = [
                        SFAPs_uV((1-all_SFAP_delays_samples(i)):end,i)
                        zeros(-all_SFAP_delays_samples(i),1)
                        ];
                end
                if (all_SFAP_delays_samples(i)>0)
                    SFAPs_uV(:,i) = [
                        zeros(all_SFAP_delays_samples(i),1)
                        SFAPs_uV(1:(end-all_SFAP_delays_samples(i)),i)
                        ];
                end
            end
            % Truncate the min_delay_samples zeros to get the intended effect of negative shifts
            SFAPs_uV = SFAPs_uV((-min_delay_samples+1):end,:);

            % Sum SFAPs to get CNAP
            CNAP_signal_uV = sum(SFAPs_uV,2);
            % save as a column vector to match CNAP signal
            final_time_vector_ms = [((1:length(CNAP_signal_uV))-1)*dt_ms_original]'; % [ms]

        end

        %{
            This function calculates the electric potentials at the specified
            coordinates using the provided z coordinates and electric
            potentials. The z coordinates must be provided as a 1-by-N vector,
            and the electric potentials must be provided as a 1-by-N vector. The
            coordinates must be at the XY origin for this particular version of
            the code, although more general versions that depend on COMSOL allow
            arbitrary coordinates.
        %}
        function [recording_sensivitity_all_fibers,n_compartments_per_fiber] = get_electric_potentials(...
                extracellular_recording_model_structure,coords_all_compartments_all_axons)

            coords_mat = cell2mat(coords_all_compartments_all_axons');


            % Create a 1-by-M vector to store forward model values
            forward_model_mat = zeros(1, size(coords_mat, 2));

            % Check whether the 3-by-K matrix consists of
            % coordinates along the same xy
            ROUND_PRECISION = 6;
            unique_xy_coords = unique(round(coords_mat([1 2],:),ROUND_PRECISION)','rows');
            if (size(unique_xy_coords,1)==1)
                one_xy = true;
                min_z = floor(min(coords_mat(3,:)));
                max_z = ceil(max(coords_mat(3,:)));
            else
                one_xy = false;
            end

            % otherwise, assert that the model structure contains a field named 'z_coords' a
            % field named 'V' so that it contains potentials along a straight line at the XY
            % origin and that all the provided coordinates satisfy this
            if (isfield(extracellular_recording_model_structure,'filename'))
                assert(strcmp(extracellular_recording_model_structure.filename(end-3:end),'.mph'),...
                    'The filename field of the extracellular recording model structure must be a COMSOL .mph file');

                import com.comsol.model.util.*
                import com.comsol.model.*
                % load the extracellular recording model
                recording_model = mphload(extracellular_recording_model_structure.filename);
                % get scaling factor relative to [m]; scale up input coords as needed
                conversion_factor = 1e-6/recording_model.geom('geom1').geomScaleFactor(); % conversion from um to COMSOL units


                % if so, set the extracellular_recording_model_structure.z_coords,...
                %                         extracellular_recording_model_structure.V
                %                         variables now; if not, then proceed to do the
                if one_xy
                    % Get the initial interpolant data on the first iteration
                    % Pass a 3-by-K matrix into the function process_matrix and store the
                    % output into the 1-by-K vector output
                    target_step_size = 1; % [um]
                    extracellular_recording_model_structure.z_coords = linspace(min_z,max_z,round((max_z-min_z)/target_step_size));
                    xyz_coords_straight_centered = [
                        unique_xy_coords(1)*ones(size(extracellular_recording_model_structure.z_coords))
                        unique_xy_coords(2)*ones(size(extracellular_recording_model_structure.z_coords))
                        extracellular_recording_model_structure.z_coords
                        ];
                    extracellular_recording_model_structure.V = mphinterp(recording_model,'V',...
                        'coord',conversion_factor*xyz_coords_straight_centered,...
                        'dataset',extracellular_recording_model_structure.dset_tag,...
                        'outersolnum',extracellular_recording_model_structure.outersolnum,'recover','pprint','matherr','on');
                end
            else
                assert(isfield(extracellular_recording_model_structure,'z_coords'),...
                    'The extracellular recording model structure must contain a field named ''z_coords''');
                assert(isfield(extracellular_recording_model_structure,'V'),...
                    'The extracellular recording model structure must contain a field named ''V''');
                assert(one_xy,'With this method, all xy coordinates must be the same such that all fibers are at the same xy location');
            end

            % Iterate through the 3-by-M matrix input in batches of 3-by-K
            % matrices, where K is batch_size; extraction occurs in batches
            % to avoid Java Out of Memory errors
            batch_size = 2e6;
            for i = 1:batch_size:size(coords_mat, 2)
                % If the remaining number of columns is less than K, then set K to the
                % remaining number of columns
                batch_size = min(batch_size, size(coords_mat, 2) - i + 1);

                coords_mat_batch_i = coords_mat(:, i:i+batch_size-1);

                % Otherwise, assert that the xyz values all lie on a straight line at the XY
                % origin, and then use interp1 to interpolate all the potentials
                if (isfield(extracellular_recording_model_structure,'z_coords'))
                    % Pass a 3-by-K matrix into the function process_matrix and store the
                    % output into the 1-by-K vector output
                    forward_model_mat_batch_i = interp1(...
                        extracellular_recording_model_structure.z_coords,...
                        extracellular_recording_model_structure.V,...
                        coords_mat_batch_i(3,:),'spline');
                else

                    % Pass a 3-by-K matrix into the function process_matrix and store the
                    % output into the 1-by-K vector output
                    forward_model_mat_batch_i = mphinterp(recording_model,'V',...
                        'coord',conversion_factor*coords_mat_batch_i,...
                        'dataset',extracellular_recording_model_structure.dset_tag,...
                        'outersolnum',extracellular_recording_model_structure.outersolnum,...
                        'recover','pprint','matherr','on');
                end
                forward_model_mat(i:i+batch_size-1) = forward_model_mat_batch_i;
            end


            % replace NaN values with zeros so that subseqeuent resampling can occur without error
            forward_model_mat(isnan(forward_model_mat)) = 0;

            % Reconvert the COMSOL output into a cell array of potentials, which is
            % more readily input to the filter calculation
            n_compartments_per_fiber = cellfun(@(x) size(x,2),coords_all_compartments_all_axons)';
            recording_sensivitity_all_fibers = mat2cell(...
                transpose(forward_model_mat),... % make each forward model a column vector
                n_compartments_per_fiber,1);


        end

        % Returns a cell array of length M containing 3-by-L_i matrices in which
        % the first, second, and third row correspond to the x, y, and z
        % coordinates (respectively) of each compartment of a given fiber among the
        % M fiber diameters. L_i is the number of compartments in total.
        % geometry_determination_method = 0 for preset fiber diameters; geometry_determination_method = 1 for MRG-based geometry interpolation; geometry_determination_method = 2 for GeometryBuilder fits from SPARC Y2Q1
        % fiber_diameters_um - a length M vector containing the fiber diameter (in
        % micrometer) of each fiber diameter
        % xy_coords - a 2-by-M matrix where the first row is the x coordinate and
        %   the second row is the y coordinate of the M axons (all coordinates in
        %   micrometer)
        %
        function [coords_all_axons,n_compartments_per_repeatable_unit,...
                distance_between_compartments_or_nodes_all_fibers] = get_straight_fiber_xyz_coords(...
                fiber_type,geometry_determination_method, fiber_diameters_um,xy_coords,LENGTH_AXONS,DESIRED_CENTER,n_stin_per_segment,flag_jitter_z)

            assert(strcmp(fiber_type,'myelinated')==1 || strcmp(fiber_type,'unmyelinated')==1, 'Invalid fiber_type value')
            % if xy_coords has only one row, assume the user specified a row vector,
            % and transpose it to make the expected column vector
            if (size(xy_coords,1)==1)
                xy_coords = xy_coords';
            end
            assert(size(xy_coords,2)==length(fiber_diameters_um) || size(xy_coords,2)==1,'you must specify one xy coordinate per fiber or only a single xy coordinate for all fibers')
            if (size(xy_coords,2)==1)
                xy_coords = repmat(xy_coords,1,length(fiber_diameters_um));
            end

            if (~exist('n_stin_per_segment','var') || isempty(n_stin_per_segment))
                n_stin_per_segment = 6; % default to six STIN per internode
            end

            if (~exist('flag_jitter_z','var') || isempty(flag_jitter_z))
                flag_jitter_z = 0;
            end

            if (strcmp(fiber_type,'myelinated'))
                n_compartments_per_repeatable_unit = n_stin_per_segment+5; % Code expects a total of one node, two MYSA, two FLUT, and however many STIN (so # STIN + 5 comparments per repeatable unit; this is asserted below
            else
                n_compartments_per_repeatable_unit = 1;
            end
            coords_all_axons = cell(length(fiber_diameters_um),1);
            distance_between_compartments_or_nodes_all_fibers = zeros(size(coords_all_axons));

            for fiber_idx = 1:length(fiber_diameters_um)

                if (strcmp(fiber_type,'myelinated'))
                    obj.fiber_dependency_method = geometry_determination_method;
                    obj.fiber_diameter = fiber_diameters_um(fiber_idx);

                    obj = CAPulator.get_fiber_geometric_props(obj);
                    stin_length = (obj.internode_length-obj.node_length-2*obj.mysa_length-2*obj.paranode_length_2)/n_stin_per_segment;

                    axon_repeatable_unit_lengths = [...
                        obj.node_length, ...
                        obj.mysa_length, ...
                        obj.paranode_length_2, ...
                        repmat(stin_length,1,n_stin_per_segment), ... % repeated by N_STIN_PER_SEGMENT
                        obj.paranode_length_2, ...
                        obj.mysa_length]; % [um]
                    INTERNODAL_DISTANCE = sum(axon_repeatable_unit_lengths);
                    assert(length(axon_repeatable_unit_lengths)==n_stin_per_segment+5,'Code expects a total of one node, two MYSA, two FLUT, and however many STIN (so # STIN + 5 comparments per repeatable unit)');
                    simulation_i_number_nodes = floor(LENGTH_AXONS/sum(axon_repeatable_unit_lengths));
                    if (mod(simulation_i_number_nodes,2)==0) % ensure number of nodes is an odd number
                        simulation_i_number_nodes = simulation_i_number_nodes +1;
                    end
                    % make the coordinates of the first repeatable piece
                    z_coords_one_axon_first_repeatable_piece = 0;
                    for i=2:length(axon_repeatable_unit_lengths)
                        z_coords_one_axon_first_repeatable_piece(i) = z_coords_one_axon_first_repeatable_piece(i-1)+...
                            axon_repeatable_unit_lengths(i-1)/2+axon_repeatable_unit_lengths(i)/2;
                    end
                    % make all the coordinates by shifting the coordinates by the internodal distance
                    % make it a column vector so that it can be added to the z_coords_one_axon_first_repeatable_piece vector
                    % via bsxfun to produce a matrix of all the z coordinates
                    shifts_to_make_all_axons = ([0:(simulation_i_number_nodes-1)])*...
                        INTERNODAL_DISTANCE;
                    z_coords_one_axon = bsxfun(@plus,shifts_to_make_all_axons',z_coords_one_axon_first_repeatable_piece);
                    % reshape the z_coords_one_axon matrix into a row vector in ascending order
                    z_coords_one_axon = reshape(transpose(z_coords_one_axon),1,[]);
                    % compartment_lengths = repmat(axon_repeatable_unit_lengths,1,simulation_i_number_nodes); % [um]
                    % cut out the trailing non-node compartments from z_coords_one_axon at the distal end so that the axon starts and ends with a nodal compartment
                    z_coords_one_axon = z_coords_one_axon(1:(end-(length(axon_repeatable_unit_lengths)-1)));
                    middle_coord_ind = (length(z_coords_one_axon)+1)/2;
                    % assert(middle_coord_ind== median(1:length(axon_repeatable_unit_lengths):length(z_coords_one_axon)),'Middle coordinate index must match old version of calculation...')
                    CENTER_NODE_Z_COORD = z_coords_one_axon(middle_coord_ind);
                    % z_coords_one_axon = z_coords_one_axon - CENTER_NODE_Z_COORD + DESIRED_CENTER;
                    z_coords_one_axon = z_coords_one_axon - CENTER_NODE_Z_COORD + DESIRED_CENTER;
                else
                    INTERNODAL_DISTANCE = 8.333; % [um]
                    simulation_i_number_nodes = round(LENGTH_AXONS/INTERNODAL_DISTANCE);
                    if (mod(simulation_i_number_nodes,2)==0) % ensure number of nodes is an odd number
                        simulation_i_number_nodes = simulation_i_number_nodes +1;
                    end
                    z_coords_one_axon = 0:INTERNODAL_DISTANCE:((simulation_i_number_nodes-1)*INTERNODAL_DISTANCE); % [um]
                    middle_coord_ind = (length(z_coords_one_axon)+1)/2;
                    CENTER_NODE_Z_COORD = z_coords_one_axon(middle_coord_ind);
                    % assert(ceil(median(1:length(z_coords_one_axon)))==middle_coord_ind)
                    z_coords_one_axon = z_coords_one_axon - CENTER_NODE_Z_COORD + DESIRED_CENTER;
                end


                %%% Calculate the distance between two adjancent compartments (unmyelinated) or
                %%% nodes (myelinated)
                distance_between_compartments_or_nodes_all_fibers(fiber_idx) = INTERNODAL_DISTANCE;

                % if specified, add jitter to z coords by up to half internodal
                % distance in either direction
                if (flag_jitter_z)
                    z_coords_one_axon = z_coords_one_axon + distance_between_compartments_or_nodes_all_fibers(fiber_idx)*(rand(1)-0.5);
                end

                % construct x & y coordinates for axon, and shift axon to center it at (0,0,0)
                coords_one_axon = zeros(3,length(z_coords_one_axon));
                coords_one_axon(1,:) = xy_coords(1,fiber_idx); %*ones(size(z_coords_one_axon));
                coords_one_axon(2,:) = xy_coords(2,fiber_idx); %*ones(size(z_coords_one_axon));
                coords_one_axon(3,:) = z_coords_one_axon;
                coords_all_axons{fiber_idx} = coords_one_axon;
            end

            for fiber_ind = 1:length(coords_all_axons)
                fiber_i_coords = coords_all_axons{fiber_ind};
                n_compartments_to_keep = length(fiber_i_coords) - mod(length(fiber_i_coords),n_compartments_per_repeatable_unit);
                coords_all_axons{fiber_ind} = fiber_i_coords(:,1:n_compartments_to_keep);
            end

        end

        %{
            This function extracts the action potential templates from the
            specified file, simplifies the data, and saves it to a new file.
        %}
        function [all_conduction_velocity_m_per_s,all_fiber_diameters,...
                common_time_vector_ms,template_data_matrix,reference_peak_times_all_templates,...
                reference_compartment_locations_um_all_templates,fiber_type] = ...
                get_template_data(template_data_source_filename,flag_dipole)
            arguments
                template_data_source_filename
                flag_dipole = 1 % use dipolar templates by default
            end

            %%% Load data
            load(template_data_source_filename,'output_data_structure','fiber_type');
            % Ensure output_data_structure is a column vector
            if (size(output_data_structure,2)>1)
                output_data_structure = transpose(output_data_structure);
            end
            % Ensure fiber_type is spelled out
            if (strcmp(fiber_type,'myel'))
                fiber_type = 'myelinated';
            end
            if (strcmp(fiber_type,'unmyel'))
                fiber_type = 'unmyelinated';
            end
            n_time_points = length(output_data_structure(1).temporal_templates);
            n_templates = length(output_data_structure);
            % identify whether to process (monopole) temporal templates or dipole
            % temporal templates
            if (flag_dipole)
                templates_name = 'dipole_temporal_templates';
                reference_peak_times_all_templates = vertcat(output_data_structure.reference_peak_time_i_dipole_ms);
            else
                templates_name = 'temporal_templates';
                reference_peak_times_all_templates = vertcat(output_data_structure.reference_peak_time_i_ms);
            end
            if (size(output_data_structure(1).(templates_name),1)==n_time_points)
                time_dim = 1;
                compartment_dim = 2;
            elseif (size(output_data_structure(1).(templates_name),2)==n_time_points)
                time_dim = 2;
                compartment_dim = 1;
            else
                error('Invalid size of output_data_structure(i).(templates_name)...');
            end

            n_compartments_per_node = size(output_data_structure(1).temporal_templates,compartment_dim);
            assert(n_templates~=n_time_points,'n_time_points and n_templates should never be equal; this is likely a bug');
            template_data_matrix = zeros(n_time_points,n_templates,n_compartments_per_node);

            for i = 1:n_templates
                % Permute so that the time dimension (the longest dimension) is
                % first
                % Make the second dimension template number
                % Make the third dimension be the MYSA, FLUT, STIN, etc. for myel (does
                % nothing for unmyel)


                template_data_matrix(:,i,:) = permute(output_data_structure(i).(templates_name),[time_dim 3 compartment_dim]);
            end


            common_time_vector_ms = output_data_structure(1).common_time_vector_ms;

            all_fiber_diameters = vertcat(output_data_structure.fiber_diameter);
            all_conduction_velocity_m_per_s = vertcat(output_data_structure.conduction_velocity_m_per_s);


            % Restore stim pulse offset to the template peak timing
            if (strcmp(fiber_type,'myelinated'))
                reference_peak_times_all_templates = reference_peak_times_all_templates + 0.1; % [ms]
            elseif (strcmp(fiber_type,'unmyelinated'))
                reference_peak_times_all_templates = reference_peak_times_all_templates + 0.4; % [ms]
            end


            % Get reference peak times and reference compartments; this will be used to
            % adjust the timing of the SFAP later on to match the source templates
            reference_compartment_locations_um_all_templates = [];
            for i = 1:size(output_data_structure,1)
                reference_compartment_locations_um_all_templates = [reference_compartment_locations_um_all_templates;...
                    1e3*output_data_structure(i).z_locations_mm_all(output_data_structure(i).target_compartment_index)];
            end

        end

        %{
            Wrapper to MATLAB's built-in fft function. It calculates a frequency vector
            based on the provided fs. MATLAB's built-in fft function is agnostic to this
            information, so it doesn't calculate freq, yet freq is very useful for
            interpreting fft results.
        %}
        function [Y,freq] = fft_wrapper(X, N, DIM, fs)

            if ~exist('DIM','var') || isempty(DIM)
                DIM = 1;
            end

            if ~exist('N','var') || isempty(N)
                N = size(X,DIM);
            end

            %%% Run FFT using MATLAB's built-in fft function and built-in inputs
            Y = fft(X,N,DIM);

            %%% Calculate freq using fs
            df = fs/N;
            if (mod(N,2)==0)
                % for even FFT, the first half starts with the 0 Hz component, while the
                % second half starts with a real (i.e., not complex-valued)  high frequency
                % component; every other value has a complex conjugate pair
                freq = [((1:(N/2))-1)*df, ((1:(N/2)))*df - (N/2 + 1)*df];

            else
                % for odd FFT, the first half starts with the 0 Hz component, and the last
                % value of the first half is the complex conjugate of the first value of
                % the second half
                freq = [((1:ceil(N/2))-1)*df, ((1:floor(N/2)))*df - ceil(N/2)*df];
            end
        end

        %{

        Description: Obtain MRG fiber geometric parameters.

        Author: Pieces were put together by Edgar Pena, using the code from Jim
        Hokanson as the basis, and using code written by Aman Aberra for method 1
        and by Nikki Pelot and Edward Liang for method 2.
        %}
        function output_obj = get_fiber_geometric_props(obj)
            %
            % Method 0: MRG present diameters
            % Method 1: Run Nikki's quadratic fits to MRG data
            % Method 2: Run Nikki & Edward's fits from SPARC project

            switch obj.fiber_dependency_method
                case 0
                    output_obj = CAPulator.helper__changePropsMethod0(obj);
                case 1
                    output_obj = CAPulator.helper__changePropsMethod1(obj);
                case 2
                    output_obj = CAPulator.helper__changePropsMethod2(obj);
                otherwise
                    error('Option #%d not recognized')
            end


        end

        function output_obj = helper__changePropsMethod0(obj)
            %
            %
            %   See MRG paper for details
            fiber_diameter_all       = [1 2 5.7      7.3     8.7     10      11.5    12.8    14      15      16];

            FIBER_INDEX = find(fiber_diameter_all == obj.fiber_diameter,1);
            if isempty(FIBER_INDEX)
                error('Unable to find specifications for given fiber size')
            end

            %ROUGH DIAMETER OUTLINE
            %--------------------------------------------------------------
            %FIBER DIAMETER > AXON DIAMETER > NODE DIAMETER
            %AXON DIAMETER = FLUT DIAMETER (PARANODE 2)
            %NODE DIAMETER = MYSA DIAMETER (PARANODE 1)


            internode_length_all     = [100 200 500      750     1000    1150    1250    1350    1400    1450    1500];
            number_lemella_all       = [15 30 80       100     110     120     130     135     140     145     150];
            %node_length             CONSTANT
            node_diameter_all        = [0.7 1.4 1.9      2.4     2.8     3.3     3.7     4.2     4.7     5.0     5.5];
            %paranode_length_1       CONSTANT
            paranode_diameter_1_all  = [0.7 1.4 1.9      2.4     2.8     3.3     3.7     4.2     4.7     5.0     5.5];
            %space_p1                CONSTANT
            paranode_length_2_all    = [5 10 35       38      40      46      50      54      56      58      60];
            paranode_diameter_2_all  = [0.8 1.6 3.4      4.6     5.8     6.9     8.1     9.2     10.4    11.5    12.7];
            %space_p2                CONSTANT
            %STIN LENGTH             DEPENDENT - delta_x_all,paranode_length_1,paranode_length_2_all,n_STIN
            axon_diameter_all        = [0.8 1.6 3.4      4.6     5.8     6.9     8.1     9.2     10.4    11.5    12.7];

            obj.internode_length     = internode_length_all(FIBER_INDEX);
            obj.number_lemella       = number_lemella_all(FIBER_INDEX);
            obj.node_diameter        = node_diameter_all(FIBER_INDEX);
            obj.paranode_diameter_1  = paranode_diameter_1_all(FIBER_INDEX);
            obj.paranode_length_2    = paranode_length_2_all(FIBER_INDEX);
            obj.paranode_diameter_2  = paranode_diameter_2_all(FIBER_INDEX);
            obj.axon_diameter        = axon_diameter_all(FIBER_INDEX);

            output_obj = obj;

        end


        function output_obj = helper__changePropsMethod1(obj)

            fd = obj.fiber_diameter;

            % Nikki fit
            % Second order polynomial using Curve Fitting Tool
            obj.number_lemella          = -0.4749*fd^2 + 16.85*fd  + -0.7648;
            obj.node_diameter           = 0.01093*fd^2 + 0.1008*fd + 1.099;
            obj.paranode_diameter_1     = obj.node_diameter;
            obj.paranode_diameter_2     = 0.02361*fd^2 + 0.3673*fd + 0.7122;
            obj.axon_diameter           = obj.paranode_diameter_2;
            obj.paranode_length_2       = -0.1652*fd^2 + 6.354*fd  + -0.2862;
            %{
            The fd cutoff for the two internodal length fits was defined by the
            intersection of the two curves.
            %}
            if (fd >= 5.643)
                obj.internode_length    = -8.215*fd^2  + 272.4*fd  + -780.2;
            else
                % Linear fit between 2 um and 5.7 um points
                obj.internode_length    =                81.08*fd  + 37.84;
            end

            output_obj = obj;

        end


        function output_obj = helper__changePropsMethod2(obj)

            fd = obj.fiber_diameter;

            obj.paranode_length_2 = -0.171 * fd^2 + 6.48 * fd + -0.935;

            fiberD_to_axonD = 4; % select fit to use for estimating axon diameter from fiber diameter
            if(fiberD_to_axonD == 0)
                %Fazan Proximal
                obj.axon_diameter = 0.553 * fd + -0.024;
            elseif(fiberD_to_axonD == 1)
                %Fazan Distal
                obj.axon_diameter = 0.688 * fd + -0.337;
            elseif(fiberD_to_axonD == 2)
                % Berthold
                obj.axon_diameter = 0.0156 * fd^2 + 0.392 * fd + 0.188;
            elseif(fiberD_to_axonD == 3)
                % Friede
                obj.axon_diameter = 0.684 * fd + 0.0821;
            elseif(fiberD_to_axonD == 4)
                %Combined Fazan Proximal, Fazan Distal, and Friede. Berthold excluded for not being vagus data
                obj.axon_diameter = 0.621 * fd - 0.121;
            end
            obj.paranode_diameter_2 = obj.axon_diameter;

            obj.node_diameter = 0.321 * obj.axon_diameter + 0.37;
            obj.paranode_diameter_1 = obj.node_diameter;

            axonD_to_number_lemella = 0; %select fit to use for estimating number of lemella from axon diameter
            if(axonD_to_number_lemella == 0)
                obj.number_lemella = floor(17.4 * obj.axon_diameter + -1.74);
            elseif (axonD_to_number_lemella == 1)
                obj.number_lemella = floor(-1.17 * obj.axon_diameter^2 + 24.9 * obj.axon_diameter + 17.7);
            end

            obj.node_length = 1;
            obj.mysa_length = 3;
            obj.internode_length = (-3.22 * fd^2 + 148 * fd + -128);

            output_obj = obj;

        end

        %{
        Return the unique values within the specified vector to the specified
        precision, and return the number of instances of each of those unique
        values at the specified precision

        Inputs
        vector_of_values - the vector of values for which you want to count the
        unique values

        desired_round_precision - (default: 6) integer indicating the precision
        (i.e., 10^-desired_round_precision, so default precision is 10^-6) with
        which the values will be considered unique; the output is rounded to the
        specified precision; negative integer values are allowed (e.g., a
        desired_round_precision value of -1 rounds numbers to the nearest 10^+1)

        Examples
        [unique_values, n_of_each_unique_value] = count_unique_values([12 9 22 22 30])
        unique_values =
            9    12    22    30
        n_of_each_unique_value =
            1
            1
            2
            1
        %}
        function [unique_values, n_of_each_unique_value,IA,IC] = ...
                count_unique_values(vector_of_values,desired_round_precision)

            if (~exist('desired_round_precision','var') || isempty(desired_round_precision))
                desired_round_precision = 6;
            end

            [unique_values,IA,IC] = unique(round(vector_of_values,desired_round_precision));
            n_of_each_unique_value = arrayfun(@(x) sum(IC==x), min(IC):max(IC))';
        end
    end
end