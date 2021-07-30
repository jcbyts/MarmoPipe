classdef NeuroNexus64_2_32_flipped < hardware.electrode.probe
    
    properties
    end
    
    methods        
        function p = NeuroNexus64_2_32_flipped(varargin)
            p@hardware.electrode.probe(varargin{:}); % base class constructor
            
            p.name = 'NeuroNexus64_2_32';
            p.manufacturer = 'NeuroNexus';
            p.design       = '64_2x32';
            p.num          = '2019?';
            p.xcoords      = [zeros(32,1) 200*ones(32,1)];
            p.ycoords      = fliplr(1:32)' * [35 35];
            p.zcoords      = zeros(32,2);
            
            
            shank1 = [1 32 2 31 3 30 4 29 ...
                5 28 6 27 7 26 8 25 ...
                9 24 10 23 11 22 12 21 ...
                13 20 14 19 15 18 16 17];
            shank2 = [33 64 34 63 35 62 36 61 ...
                37 60 38 59 39 58 40 57 ...
                41 56 42 55 43 54 44 53 ...
                45 52 46 51 47 50 48 49];


            % Combine shank1 and 2 into one channel map
            shank = [shank1 shank2];
            
            % full channel map for the omnetics36 connector (36 channels)
            % use NAN for the Ground and Reference Channels
            omnetics1 = [nan 37 39 40 42 43 45 46 48 17 19 20 22 23 25 26 28 nan ...
                nan 36 38 35 41 34 44 33 47 18 32 21 31 24 30 27 29 nan];
            omnetics2 = [nan 64 62 60 58 56 54 52 50 15 13 11 9 7 5 3 1 nan ...
                nan 63 61 59 57 55 53 51 49 16 14 12 10 8 6 4 2 nan];
            
            % combine both omnetics connectors (72 channels) Flip one
            % connector
            omnetics = [omnetics1 rot90(omnetics2,2)];
            
            % reshape to match the physical shape of the connectors
            omnetics = reshape(omnetics', [], 4);

            % loop over channels and find the map from probe to omnetics
            % connector
            Nchan = numel(shank);
            p.channelMap = nan(1,Nchan);
            for i = 1:Nchan
                p.channelMap(i) = find(omnetics==shank(i));
            end
            
            p.connector    = 'Omnetics36';
            p.material     = '??';
            p.impedence    = 0.4;
            
        end
        
    end
end