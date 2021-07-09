function [tstrobes, strobes] = read_vpx_strobes( filename )
% function [tstrobes, strobes] = read_vpx_strobes( filename )
%***
%*** takes a VPX file, searches for lines tagged by 12, and
%*** identifies TRIALSTART or TRIALEND as the tag
%*** then reads the six numbered taglet on that line
%*** encodes everything like ephys strobes and tstrobes so 
%*** tstrobes are timestamps in vpx time
%*** and strobes are (63 start, 62 end) followed by tags in a Nx1 column
%***
%*** input:   filename
%***
%*** output:  tstrobes and strobes
%***       - tstrobes are timestamps in vpx time
%***       - strobes are integers from 0 to 63 (like ephys strobes)
%************************************************

  %***** read in the raw eye data first
  MaxCnt = 5000000;  % max, 5 million samples, thats 1.5 hrs at 1000 Hz
  tstrobes = zeros(MaxCnt,1);
  strobes = uint8(zeros(MaxCnt,1));
  START_TAG = 63; %64;
  END_TAG = 62; %63;
  STARTPAT = 'TRIALSTART:TRIALNO:';
  STARTN = size(STARTPAT,2);
  ENDPAT = 'TRIALENDED:TRIALNO:';
  ENDN = size(ENDPAT,2);
  %*******************
  disp(sprintf('Reading start and end trial strobes from VPX %s',filename));
  %*************
  fd = fopen(filename,'r');
  if (fd)
    tline = fgets(fd);
    cnt = 0;
    while (tline ~= -1)
       ltag = str2num(tline(1:2));
       if (ltag == 12) %strobe tag, this is eye position
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
              tt = dval(2); %timestamp
              cnt = cnt + 1;
              tstrobes(cnt) = tt;
              strobes(cnt) = uint8(tag);
              for k = 1:6   % write the taglet
                  tstrobes(cnt+k) = tt;
                  strobes(cnt+k) = uint8( dval(3+k) );
              end
              cnt = cnt + 6;
              %********* update display
              if (mod(cnt,140)==0)
                  disp(sprintf('VPX Strobe Count %d',cnt));
              end
           end
           %****************
           
       end
       %*************
       if (cnt > MaxCnt)
           disp('WARNING: Truncating data read, reached max count');
           break;
       end
       %***********
       tline = fgets(fd);  % read the next line
       %***********
    end
    fclose(fd);
  end  % if file
  %******* when all done, throw away extra zeros in data structs
  if (cnt)
    tstrobes = tstrobes(1:cnt,:);
    strobes = strobes(1:cnt,:);
  else
    disp('Warning: no strobes found in VPX file');
    tstrobes = [];
    strobes = [];
  end
  %*****************

return;