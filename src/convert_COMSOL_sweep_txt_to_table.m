%{

Author: Edgar Peña

Given a text file formatted to be loaded into COMSOL for a parameteric or batch
sweep, create a MATLAB table from the data in the text file. The data in COMSOL
parameteric sweep or batch sweep files have the following format:

    parameter1 "value1, value2, value3, ..., valueN" [unit_of_parameter1]
    parameter2 "value1, value2, value3, ..., valueN" [unit_of_parameter2]
    ...
    parameterN "value1, value2, value3, ..., valueN" [unit_of_parameterN]

%}
function [parameter_values_table,parameter_units_table] = ...
    convert_COMSOL_sweep_txt_to_table(filename)
    % read the file
    file = fileread(filename);
    
    % split the file into lines
    lines = splitlines(file);
    
    % create a cell array to store the data
    data_struct = [];
    
    % iterate through each line
    for i = 1:length(lines)
        
        % skip empty line at the end of the file
        if (isempty(lines{i}) && i==length(lines))
            % do nothing

        else
            % get rid of all spaces, since they end up in the parameter
            % name which produces invalid structure fields
            lines{i} = strrep(lines{i},' ','');

            % split the line into the parameter and the values
            splitLine = split(lines{i}, '"');

            % get the parameter name
            parameter = splitLine{1};

            % get the values
            values = split(splitLine{2}, ',');
            % extract numeric values from text, then convert to vector
            % of doubles rather than cell array
            for value_ind = 1:length(values)
                values{value_ind} = sscanf(values{value_ind},'%f');
            end
            values = cell2mat(values);


            % get the units
            units = split(splitLine{3}, '[');
            units = units{2};
            units = split(units, ']');
            units = units{1};

            % store the data in the cell array
            for value_ind = 1:length(values)
                data_struct(value_ind).(parameter) = values(value_ind);
            end
            units_struct(i).ParameterName = parameter;
            units_struct(i).ParameterUnits = units;
        end
    end
    
    % Create a table from the data structure such that there is one column for each
    % parameter, and one row for each value. 
    parameter_values_table = struct2table(data_struct);

    % Also create a separate table from the units structure that maps the
    % parameter names to the parameter units such that there are two
    % columsn in the table: ParameterName, ParameterUnits
    parameter_units_table = struct2table(units_struct);

end