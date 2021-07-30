function [data, meta_file] = getMetaTable(meta_file)
% % Get Meta Table reads in the experimental meta data and returns a table
% meta = getMetaTable()

if nargin < 1
    % SERVER_DATA_DIR = getpref('EPHYS', 'SERVER_DATA');
    SERVER_DATA_DIR = 'Z:\Data\PLDAPS\Ellie';

    meta_file = fullfile(SERVER_DATA_DIR, 'meta_data.xls');
end

[~, ~, alldata] = xlsread(meta_file);
    
data = cell2table(alldata(2:end,:), 'VariableNames', alldata(1,:));