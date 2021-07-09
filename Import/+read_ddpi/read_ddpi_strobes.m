function [tstrobes, strobes] = read_ddpi_strobes( filename )
% function [tstrobes, strobes] = read_ddpi_strobes( filename )
%***
%*** takes a DDPI file, searches for lines with messages non-zero
%*** and decodes their time tags, reverting them back to six numbered taglets
%*** which can be matched to the Matlab start and stop clock tags
%***
%*** encodes everything like ephys strobes and tstrobes so 
%*** tstrobes are timestamps in vpx time
%*** and strobes are (63 start, 62 end) followed by tags in a Nx1 column
%***
%*** input:   filename
%***
%*** output:  tstrobes and strobes
%***         - tstrobes are timestamps in vpx time
%***         - strobes are integers from 0 to 63 (like ephys strobes)
%************************************************

  %***** read in the raw eye data first
  Output = read_ddpi.ddpiReadFile(filename);
  % 1) signalType;
  % 2) time;
  % 3) p1x;
  % 4) p1y;
  % 5) p1r;
  % 6) p1I;
  % 7) p4x;
  % 8) p4y;
  % 9) p4r;
  % 10) p4I;
  % 11) p4score;
  % 12) tag;
  % 13) message;
  %***********
  
  %******* First identify all the clock messages (will be huge number
  % because Jake takes the 6 CLOCK integers and combines to make a
  % 12 digit integer from them ********************
  zclock = find( Output(13,:) > 10e6 );
  tstrobes = [];
  strobes = [];
  if isempty(zclock)
      disp('Warning: no strobes found in DDPI file');
      return;
  end
  
  %******** OK, the next part is unfortunately a bit complicated, but
  %******** the odd design comes from the DataPixx trial coding, that's why
  %****** you make a list of strobes and their times:
  %****** the first strobe is a start or end code (START_TAG or END_TAG)
  %******  then comes six numbers, which are the taglet
  %****** and ideally the END_TAG will follow a START_TAG
  
  % now search those codes to categorize as start or end
  disp(sprintf('Reading start and end trial strobes from DDPI %s',filename));
  %************
  START_TAG = 63; %64;
  END_TAG = 62; %63;
  for k = 1:length(zclock)
     tt = Output(2,zclock(k)) / 1000;  %convert ms to secs (consistent VPX)
     if (zclock(k)+20) > size(Output,2)
         continue;
     end
     postval = mean( Output(1,(zclock(k)+1):(zclock(k)+20)) );  % tag is non-zero in trial
     if (postval == 0)
        tag = END_TAG; 
     else
        tag = START_TAG;
     end
     %****** add the code back to the list
     tstrobes = [tstrobes ; tt];  % time of code
     strobes = [strobes ; tag];
     %*** make a taglet list, in reverse and swap back after
     codon = Output(13,zclock(k));
     tlist = [];
     for j = 1:6   % write out the taglet from 12 digit code
          val = mod(codon,100);
          tlist = [tlist ; val];
          codon = floor( codon/100 );
     end
     tlist = flipud(tlist);
     for j = 1:6
          tstrobes = [tstrobes ; tt];
          strobes = [strobes ; tlist(j)];
     end
     %*****************
  end
  strobes = uint8(strobes);
  disp(sprintf('Completed reading strobes from DDPI %s',filename));
  
  
  %***** read in the raw eye data first
%   MaxCnt = 5000000;  % max, 5 million samples, thats 1.5 hrs at 1000 Hz
%   tstrobes = zeros(MaxCnt,1);
%   strobes = uint8(zeros(MaxCnt,1));
%   START_TAG = 63; %64;
%   END_TAG = 62; %63;
%   STARTPAT = 'TRIALSTART:TRIALNO:';
%   STARTN = size(STARTPAT,2);
%   ENDPAT = 'TRIALENDED:TRIALNO:';
%   ENDN = size(ENDPAT,2);
%   %*******************
%   disp(sprintf('Reading start and end trial strobes from VPX %s',filename));
%   %*************
%   fd = fopen(filename,'r');
%   if (fd)
%     tline = fgets(fd);
%     cnt = 0;
%     while (tline ~= -1)
%        ltag = str2num(tline(1:2));
%        if (ltag == 12) %strobe tag, this is eye position
%            kstart = strfind(tline,STARTPAT);
%            kend = strfind(tline,ENDPAT);
%            tag = 0;  % default skip the line if not identified as start or end
%            if ~isempty(kstart)
%               tline2 = tline([1:(kstart-1),(kstart+STARTN):end]); 
%               tag = START_TAG;
%            end
%            if ~isempty(kend)
%               tline2 = tline([1:(kend-1),(kend+ENDN):end]); 
%               tag = END_TAG; 
%            end
%            %******** if identified, then write out strobes
%            if (tag > 0)
%               dval = sscanf(tline2,'%f');
%               tt = dval(2); %timestamp
%               cnt = cnt + 1;
%               tstrobes(cnt) = tt;
%               strobes(cnt) = uint8(tag);
%               for k = 1:6   % write the taglet
%                   tstrobes(cnt+k) = tt;
%                   strobes(cnt+k) = uint8( dval(3+k) );
%               end
%               cnt = cnt + 6;
%               %********* update display
%               if (mod(cnt,140)==0)
%                   disp(sprintf('VPX Strobe Count %d',cnt));
%               end
%            end
%            %****************
%            
%        end
%        %*************
%        if (cnt > MaxCnt)
%            disp('WARNING: Truncating data read, reached max count');
%            break;
%        end
%        %***********
%        tline = fgets(fd);  % read the next line
%        %***********
%     end
%     fclose(fd);
%   end  % if file
%   %******* when all done, throw away extra zeros in data structs
%   if (cnt)
%     tstrobes = tstrobes(1:cnt,:);
%     strobes = strobes(1:cnt,:);
%   else
%     disp('Warning: no strobes found in VPX file');
%     tstrobes = [];
%     strobes = [];
%   end
%   %*****************

return;