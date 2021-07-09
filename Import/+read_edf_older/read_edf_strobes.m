function [tstrobes, strobes] = read_edf_strobes( FEVENTS )
% function [tstrobes, strobes] = read_edf_strobes( FEVENTS )
%***
%*** takes a EDF matlab FEVENTS Struct, lines events and
%*** identifies TRIALSTART or TRIALEND as the tag
%*** then reads the six numbered taglet on that line
%*** encodes everything like ephys strobes and tstrobes so 
%*** tstrobes are timestamps in vpx time
%*** and strobes are (63 start, 62 end) followed by tags in a Nx1 column
%***
%*** input:   FEVENTS struct, from the Matlab EDF struct
%***
%*** output:  tstrobes and strobes
%***       - tstrobes are timestamps in vpx time
%***       - strobes are integers from 0 to 63 (like ephys strobes)
%************************************************

  %***** read in the raw eye data first
  MaxCnt = 5000000;  % max, 5 million samples, thats 1.5 hrs at 1000 Hz
  tstrobes = zeros(MaxCnt,1);
  strobes = uint8(zeros(MaxCnt,1));
  START_TAG = 63; 
  END_TAG = 62; 
  STARTPAT = 'TRIALSTART:TRIALNO:';
  STARTN = size(STARTPAT,2);
  ENDPAT = 'TRIALENDED:TRIALNO:';
  ENDN = size(ENDPAT,2);
  %*******************
  disp(sprintf('Reading start and end trial strobes from EDF struct '));
  %*************
  cnt = 0;
  for i = 1:size(FEVENTS,2)
    if (FEVENTS(i).type == 24) % message event
        tline = FEVENTS(i).message;  
        kstart = strfind(tline,STARTPAT);
        kend = strfind(tline,ENDPAT);
        tag = 0;  % default skip the line if not identified as start or end
        if ~isempty(kstart)
           tline2 = tline([1:(kstart-1),(kstart+STARTN):end]); 
           tag = START_TAG;
        end
        if ~isempty(kend)
           tline2 = tline([1:(kend-1),(kend+ENDN):end]); 
           tag = END_TAG; 
        end
        %******** if identified, then write out strobes
        if (tag > 0)
              dval = sscanf(tline2,'%f');
              tt = (double(FEVENTS(i).sttime)/1000); % dval(2); %timestamp
              cnt = cnt + 1;
              tstrobes(cnt) = tt;
              strobes(cnt) = uint8(tag);
              for k = 1:6   % write the taglet
                  tstrobes(cnt+k) = tt;
                  strobes(cnt+k) = uint8( dval(1+k) );  %dval(3+k)
              end
              cnt = cnt + 6;
              %********* update display
              if (mod(cnt,120)==0)
                  disp(sprintf('EDF Strobe Count %d',cnt));
              end
              %**************
        end
        %****************    
        if (cnt > MaxCnt)
           disp('WARNING: Truncating data read, reached max count');
           break;
        end
        %*********
    end
  end
  %******* when all done, throw away extra zeros in data structs
  if (cnt)
    tstrobes = tstrobes(1:cnt,:);
    strobes = strobes(1:cnt,:);
  else
    disp('Warning: no strobes found in EDF file');
    tstrobes = [];
    strobes = [];
  end
  %*****************

return;
