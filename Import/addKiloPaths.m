function varargout = addKiloPaths(user)
% set paths for RIGTESTING project
% this assumes you are running from the RIGTESTING folder
if nargin < 1
    user = [];
end

switch user
    case 'gravedigger'
        
        % we need the full marmopipe / import paths
        marmoPipePath = 'C:\Users\Jake\Dropbox\Marmo Lab Website\PSA\Code';
        
        
        % set preferences for where the data live. Storing these paths as
        % preferences means we can recover them from any other function.
        % It's super useful.
        
        % Raw unprocessed data:
        setpref('KILOSORT', 'SERVER_DATA_DIR', 'C:\Raw')
        
        % where the data live
        dataPath = getpref('KILOSORT', 'SERVER_DATA_DIR');
        
        % processed data:
        setpref('KILOSORT', 'PROCESSED_DATA_DIR', 'C:\Processed')
        
    otherwise
        error('This must be run on a machine that has set up Kilosort 2')
end


addpath(marmoPipePath)
addMarmoPipe

if nargout == 1
    varargout{1} = dataPath;
end
    