function [eventChannels, eventSamples, evinfo] = load_kwe(filename, varargin)
% helper function for reading KWIK format event files

RecNumber = 0; % hard coded -- I think this is right

recordings = h5read(filename, '/event_types/TTL/events/recording');
timestamps = h5read(filename, '/event_types/TTL/events/time_samples');
            
eventChannels = h5read(filename, '/event_types/TTL/events/user_data/event_channels');
eventID = h5read(filename, '/event_types/TTL/events/user_data/eventID');
nodeID = h5read(filename, '/event_types/TTL/events/user_data/nodeID');

evinfo = struct('eventId', eventID, 'recordings', recordings, 'nodeID', nodeID, 'timestamps', timestamps);

%eventSamples depends on recording number.
location='/recordings';
info=h5info(filename,location);
nRecs=length(info.Groups);
for iRec=1:nRecs
    iRecNumber=(info.Groups(iRec).Name(13:end));
    if str2double(iRecNumber)==RecNumber
        %toutDir=([outDir '_R' iRecNumber]);
        
        st_index=strcmp('start_time',{info.Groups(iRec).Attributes.Name});
        recStartTime = double(info.Groups(iRec).Attributes(st_index).Value);
        eventSamples=timestamps-recStartTime;
    end
end