%------------------------------------------------------------------------------
%	Date 		: Oct 25, 2017
%	Description :
%		This file is used to make the data set for comparision algorithm to
%	recognize time and location from a string. Then save the data set to a file.
%------------------------------------------------------------------------------


%% Clean
clear
clc


%% Dataset
% Time, date
dtset.time{1,1} = 'jan ';
dtset.time{2,1} = 'feb ';
dtset.time{3,1} = 'mar ';
dtset.time{4,1} = 'apr ';
dtset.time{5,1} = 'may ';
dtset.time{6,1} = 'jun ';
dtset.time{7,1} = 'jul ';
dtset.time{8,1} = 'aug ';
dtset.time{9,1} = 'sep ';
dtset.time{10,1} = 'oct ';
dtset.time{11,1} = 'nov ';
dtset.time{12,1} = 'dec ';

dtset.time{13,1} = 'january ';
dtset.time{14,1} = 'february ';
dtset.time{15,1} = 'march ';
dtset.time{16,1} = 'april ';
dtset.time{17,1} = 'may ';
dtset.time{18,1} = 'june ';
dtset.time{19,1} = 'july ';
dtset.time{20,1} = 'august ';
dtset.time{21,1} = 'september ';
dtset.time{22,1} = 'october ';
dtset.time{23,1} = 'november ';
dtset.time{24,1} = 'december ';

% Time, hour
dtset.hour{1,1} = 'pm';
dtset.hour{2,1} = 'p.m';
dtset.hour{3,1} = 'am';
dtset.hour{4,1} = 'a.m';

% Location
dtset.location{1,1} = 'street';
dtset.location{2,1} = 'st.';
dtset.location{3,1} = 'road';
dtset.location{4,1} = 'room';
dtset.location{5,1} = 'court';
dtset.location{6,1} = 'parkway';
dtset.location{7,1} = 'pike';
dtset.location{8,1} = 'university';
dtset.location{9,1} = 'mountain';
dtset.location{9,1} = 'st,';


%% Save
save dt_set.mat

