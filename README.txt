This repository contains all the code necessary to generate and plot the figures in, "Computational Models of Compound Nerve Action Potentials: Efficient Filter-Based Methods to Quantify Effects of Tissue Conductivities, Conduction Distance, and Nerve Fiber Parameters", *except* for the following data files, which are uploaded to a data repository instead (available at https://doi.org/10.7924/r4pc3624h ):
	files that belong in bin/ folder
		bin/template_data_20221004_0_2022100409.mat
		bin/template_data_20221004_0_2022100410.mat
	files that belong in results/ folder
		results/in_vivo_data_table.mat
		results/myelinated_fiber_CNAPs_across_materials.mat
		results/unmyelinated_fiber_CNAPs_across_materials.mat


This repository does not provide code to run NEURON or to generate action potential templates (e.g., to reproduce the brute force analysis); instead, it is necessary to use the ASCENT software version released with the published manuscript to simulate nerve recordings and to generate action potential templates. Access to a computing cluster and a COMSOL license are also required. Instructions for using this functionality of ASCENT are available on the ASCENT webpage.

The list below describes each of the key code files. Data files include in vivo maximal compound nerve action potential signals (CNAPs) and associated metadata from rat cervical vagus nerve, as well as computationally generated CNAPs. For details on how the data were generated, please see the associated manuscript.

To run the code, set the current working directory to be src, then run RUN_simulate_CAP.m. That is the master script that runs everything, and the entire analysis that it conducts can be run on a normal desktop computer with a MATLAB license (R2019b or later). RUN_simulate_CAP.m calls the methods of CAPulator.m, which in turn takes its inputs from the JSON files present in code/JSON_input_params. Changing the parameters in those files enables changing the simulation (e.g., changing the fiber diameters present, the conduction distance, etc.). 

Below is a description of the key files for conducting the analysis in the manuscript:
- src/RUN_simulate_CAP.m: This generates the data for and plots every figure in the main manuscript text. It simulates CAP signals across different conduction distances and fiber parameters. The only exception is that it does not run the simulations across tissue conductivity values; instead, the output of those is loaded here and plotted. To run a subset of the simulations/data, edit the zeros and ones in the 'function_run_status_and_handle' variable.

- src/CAPulator.m: This is the main calculation code that calculates the CNAP signals for the specified params. The following methods are the key methods are contained in this class:
		calculate_CAP_via_temporal_templates - calculates the CNAP signal due to the specified parameters, where parameters are specified by the files in the JSON_input_params folder.
		get_template_data - gets the template data from the file specified in the corresponding JSON_input_params file
		get_straight_fiber_xyz_coords - generates xyz coordinates of nerve fibers based on parameters specified in the corresponding JSON_input_params file
		get_electric_potentials - gets the electric potentials at the specified locations due to the extracellular recording model specified in the corresponding JSON_input_params file

Results from our brute force analysis are saved in the following files:
- ../results/ground_truth_CNAP_myel.mat
- ../results/ground_truth_CNAP_unmyel.mat

The electric potentials for the recording sensitivity maps were extracted via COMSOL. While a more general version of src/CAPulator.m is able to use COMSOL with Livelink for MATLAB to extract potentials from COMSOL files, the version of src/CAPulator.m shared for review does not require a COMSOL license or MATLAB with LiveLink. Instead, the potentials are extracted for and stored in the following files:
- ../bin/default_volume_conductor_potentials.mat
- ../bin/tuned_volume_conductor_potentials.mat
- ../bin/tuned_volume_conductor_potentials_ch2.mat
- ../bin/tuned_volume_conductor_potentials_ch3.mat

JSON inputs to simulate CNAPS for a baseline comparison and after tuning are below. The src/RUN_simulate_CAP.m script modifies the parameters after loading them to conduct the various senstivity analyses.
    JSON_input_params/RUN_CNAP_20221004_ASCENT_myel.json 
    JSON_input_params/RUN_CNAP_20221004_ASCENT_unmyel.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_myel_11mm_ch1_polarity1.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_myel_11mm_ch2_polarity1.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_myel_11mm_ch3_polarity1.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_11mm_ch1_polarity1.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_11mm_ch2_polarity1.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_11mm_ch3_polarity1.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_myel_15mm_ch1_polarity1.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_15mm_ch1_polarity1.json
    JSON_input_params/CNAP_20221004_cylindrical_cuff_unmyel_6mm_ch1_polarity1.json
	
Other supporting files include the following:
	- inverse_transform_sample_histogram.m
	- process_params.m
	- convert_COMSOL_sweep_txt_to_table.m
	- ../results/sens_analysis_full.txt 

