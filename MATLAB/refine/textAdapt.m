%------------------------------------------------------------------------------
%	Date 		: Oct 26, 2017
%	Description :
%		This file is used to filter the sub-string that may contain the key
%	information from a string.
%------------------------------------------------------------------------------
function str = textAdapt(fn_txtOCRfinal, dataset)

%% Open file
fid = fopen(fn_txtOCRfinal);
tline = fgetl(fid);

%% Get fields of struct
names = fieldnames(dataset);


%% Scan file
flgDate = false;
str = [];
while ischar(tline)
    % Convert to lowercase
    str_lower = lower(tline);
    % Compare to dataset
    if(length(str_lower) < 60)
    for j = 1:length(names)
        data_set = getfield(dataset, names{j});
        for i = 1 : length(data_set)
            ind = strfind(str_lower, data_set{i});
            if ~isempty(ind)
                flgDate = true;
                break;
            end
        end
        if flgDate
            break;
        end
    end
    end
    % If is date
    if flgDate
        str = sprintf('%s %s\n', str, tline);
        flgDate = false;
    end
    % Read new line
    tline = fgetl(fid);
end


%% Finally, replace special substring
str = strrep(str, '—', ' - ');


%% Close file
fclose(fid);


end

