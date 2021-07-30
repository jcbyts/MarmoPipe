classdef NeuroNexus64_2_32 < hardware.electrode.probe
    
    properties
    end
    
    methods        
        function p = NeuroNexus64_2_32(varargin)
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
            shank2 = shank1 + 32;
%             shank2 = [33 64 34 63 35 62 36 61 ...
%                 37 60 38 59 39 58 40 57 ...
%                 41 56 42 55 43 54 44 53 ...
%                 45 52 46 51 47 50 48 49];

            % Combine shank1 and 2 into one channel map
            shank = [shank1 shank2];
            
            % full channel map for the omnetics36 connector (36 channels)
            % use NAN for the Ground and Reference Channels
            omnetics1 = [nan 37 39 40 42 43 45 46 48 17 19 20 22 23 25 26 28 nan ...
                nan 36 38 35 41 34 44 33 47 18 32 21 31 24 30 27 29 nan];
%             omnetics2 = [nan 64 62 60 58 56 54 52 50 15 13 11 9 7 5 3 1 nan ...
%                 nan 63 61 59 57 55 53 51 49 16 14 12 10 8 6 4 2 nan];
%             omnetics1 = [nan 34 43:45 33 46:48 17:19 32 20:22 31 nan ...
%                 nan 42:-1:40 35 39:-1:36 29:-1:26 30 25:-1:23 nan];
            omnetics2 = [nan 64:-2:50 15:-2:1 nan ...
                nan 63:-2:49 16:-2:2 nan];
%             omnetics1 = [nan 47:-2:33 32:-2:18 nan ...
%                 nan 48:-2:34 31:-2:17 nan];
         
%         
            
           % combine both omnetics connectors (72 channels)
            omnetics = [omnetics1 omnetics2];
            
            % reshape to match the physical shape of the connectors
            omnetics = reshape(omnetics', [], 4);

            % loop over channels and find the map from probe to omnetics
            % connector
            Nchan = numel(shank);
            p.channelMap = nan(1,Nchan);
            for i = 1:Nchan
                p.channelMap(i) = find(omnetics==shank(i));
            end
            
%             % the probe is 32 x 2 channels while each omnetics connector is
%             % 36 channels. There are 4 extra channels for ground and
%             % reference. We have to account for them here.
%             probe2omnetics = find(~isnan(omnetics));
%             
%             % reduce the omnetics channel map to 32 channels per connector
%             omnetics = omnetics(probe2omnetics);
%             
%             % Finally, combine all those steps
%             p.channelMap   = probe2omnetics(omnetics(shank));
            
            p.connector    = 'Omnetics36';
            p.material     = '??';
            p.impedence    = 0.4;
            
        end
        
    end
end