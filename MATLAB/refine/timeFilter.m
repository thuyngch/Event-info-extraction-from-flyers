%------------------------------------------------------------------------------
%	Date 		: Oct 28, 2017
%	Description :
%		This file is used to extract time (date and hour) sub-string from a string.
%------------------------------------------------------------------------------
function str_time = timeFilter(str_raw, dt_set)

%% Split the string
id_newline = strfind(str_raw, sprintf('\n'));
num_line = length(id_newline);
if num_line < 1
    str_time = '';
    return;
end
str = cell(num_line, 1);
str{1,1} = str_raw(1: id_newline(1)-1);
for i = 2: num_line
    str{i,1} = str_raw(id_newline(i-1)+1 : id_newline(i)-1);
end


%% Compare to the dataset
len = length(dt_set);
str_time = '';
flgTime = false;
for line = 1:num_line
    % Scan
    for i = 1:len
        % Convert to lowercase
        str_lower = lower(str{line,1});

        % Compare to the dataset
        ind = strfind(str_lower, dt_set{i});
        if ~isempty(ind)
            flgTime = true;
            break;
        end

        % Post filter
        ind = strfind(str_lower, ':');
        if ~isempty(ind)
        for j = 1:length(ind)
            low = ind(j)-1;  if(low<1); break; end
            up = ind(j)+1; if(up>length(str_lower)); break; end
        if ~isempty(str2num(str_lower(low))) && ~isempty(str2num(str_lower(up)))
            flgTime = true;
            break;
        end
        end
        end
    end

    % If valid
    if flgTime
        str_time = sprintf('%s %s\n', str_time, str{line,1});
        flgTime = false;
    end
end


end