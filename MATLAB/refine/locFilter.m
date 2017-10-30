%------------------------------------------------------------------------------
%	Date 		: Oct 28, 2017
%	Description :
%		This file is used to extract location sub-string from a string.
%------------------------------------------------------------------------------
function str_loc = locFilter(str_raw, dt_set)

%% Split the string
id_newline = strfind(str_raw, sprintf('\n'));
num_line = length(id_newline);
if num_line < 1
    str_loc = '';
    return;
end
str = cell(num_line, 1);
str{1,1} = str_raw(1: id_newline(1)-1);
for i = 2: num_line
    str{i,1} = str_raw(id_newline(i-1)+1 : id_newline(i)-1);
end


%% Compare to the dataset
len = length(dt_set);
str_loc = '';
flgLoc = false;
for line = 1:num_line
    % Scan
    for i = 1:len
        % Convert to lowercase
        str_lower = lower(str{line,1});

        % Compare to the dataset
        ind = strfind(str_lower, dt_set{i});
        if ~isempty(ind)
            flgLoc = true;
            break;
        end
    end

    % If valid
    if flgLoc
        str_loc = sprintf('%s %s\n', str_loc, str{line,1});
        flgLoc = false;
    end
end


end