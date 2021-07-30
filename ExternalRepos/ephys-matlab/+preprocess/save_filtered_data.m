function save_filtered_data(ops, overwrite)
% fproc

if nargin < 2
    overwrite = false;
end

if ~overwrite && exist(ops.fproc, 'file')==2
    disp('file already exists')
    return
end

fid = fopen(ops.fbinary);
fHP = strrep(ops.fbinary, '.dat', '_HP.dat');
if exist(fHP, 'file')==2
    delete(fHP)
end
fout = fopen(fHP, 'w');

fseek(fid, 0, 'eof');
filesize = ftell(fid);

fseek(fid, 0, 'bof');

Nchan = ops.Nchan;

nTotSamp = filesize/Nchan;

nBatch = 10;
batchSize = floor(nTotSamp/nBatch);

% highpass filter
Fs = ops.fs;
[b,a] = butter(5, (ops.fshigh/Fs)*2, 'high');

for iBatch = 1:nBatch
    if iBatch == nBatch
        bufferSize = [Nchan batchSize + mod(nTotSamp, batchSize)];
    else
        bufferSize = [Nchan batchSize];
    end
    
    % read data
    data = double(fread(fid, bufferSize, '*int16'));
    
    if isfield(ops, 'deadChannels')
        goodChannels = setdiff(1:ops.Nchan, ops.deadChannels);
        
        if size(data,2)<=1
            
        else
            data = data - mean(data(goodChannels,:));
        end
    else
        data = data - mean(data,1);
    end
    
    data = filter(b,a,data');
%     data = gpufilter(data, ops, 1:Nchan);
    
    fwrite(fout, int16((data))', 'int16');
   
end

% close up shop
fclose(fid);
fclose(fout);

fprintf('Done\n')