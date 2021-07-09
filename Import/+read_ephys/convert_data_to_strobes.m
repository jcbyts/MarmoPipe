function [tvalue,value] = covert_data_to_strobes(data,ts,info)
%****** takes inputs from the read .events file
%****** returns the strobes assuming PLDAPS strobe function format
   bitNumber=[];
   timestamps=[];
   highlow=[];

   %****** coversion coming from sync2OeClock code
   tmp_bitNumber = data;
   tmp_timestamps = ts;
   tmp_highlow=info.eventId;
   timestamps=[timestamps; tmp_timestamps(:)]; %#ok<*AGROW>
   bitNumber=[bitNumber; tmp_bitNumber(:)];
   highlow=[highlow; tmp_highlow(:)];

   strobeSet=find(bitNumber==7 & highlow==1);
   strobeUnset=find(bitNumber==7 & highlow==0);
   strobeUnset=[1; strobeUnset];
   % extract strobe values
   value=nan(size(strobeSet));
   tvalue = nan(size(strobeSet));
   for iStrobe=1:length(strobeSet)
     ts=timestamps <= timestamps(strobeSet(iStrobe)) & timestamps >= timestamps(strobeUnset(iStrobe)) & bitNumber~=7;
     value(iStrobe)=sum(2.^bitNumber(ts) .* highlow(ts));
     tvalue(iStrobe) = timestamps(strobeSet(iStrobe));
   end

return;

