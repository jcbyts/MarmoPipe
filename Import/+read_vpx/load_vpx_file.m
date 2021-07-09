function vpx  = load_vpx_file( filename )
% function  vpx = load_vpx_file( filename )
% - reads vpx file out for raw data, median filters and smooths
%   to form a file up-sampled to 1000 hz
% - also detects trial start and end strobe times with taglets
%***
%  input:  filename
%***
%  outputs:   vpx struct
%             vpx.raw - [Nx4] - fields: time, x, y, pupil
%             vpx.smo - [Tx4] - same fields, but upsampled 1000hz
%                               via linear interp, and smoothing with
%                               a Gaussian of sigma 5ms
%             returns [] if file problem
%
%***** it reads the VPX file strobes and reformats them
%***** to be just like ephys strobes, so we can treat them the same
%             vpx.tstrobes - timestamps of strobes
%             vpx.strobes - strobe values
%*******************************************

  MEDFILT = 2;   % median filtering of data series, 0 if none
  GSIG = 5;      % Gaussian smoothing in ms (must be odd integer), 0 if none
  MaxCnt = 5000000;  % max, 5 million samples, thats 1.5 hrs at 1000 Hz
  
  %*******************
  disp(sprintf('Reading VPX file %s',filename));
  vpx = struct;
  vpx_strobes = struct;
  
  %***** read in the raw eye data first
  dat = zeros(MaxCnt,4);  %allocate space to load
  fd = fopen(filename,'r');
  if (fd)
   tline = fgets(fd);
   cnt = 1;
   while (tline ~= -1)
       ltag = str2num(tline(1:2));
       if (ltag == 10) %datadum, this is eye position
           dval = sscanf(tline,'%f');    
           if (size(dval,1) >= 5)
               xx = dval(4);  %(dval(4)-Ac1)/(Adx*SpixPerDeg);
               yy = dval(5);  %((1-dval(5))-Ac2)/(Ady*SpixPerDeg);
               tt = dval(2);
               pp = (dval(9)+dval(10))/2;
               dat(cnt,:) = [tt xx yy pp];  
               cnt = cnt + 1;
           end
       end
       %*************
       if (cnt > MaxCnt)
           disp('WARNING: Truncating data read, reached max count');
           break;
       end
       %**********
       if (mod(cnt,50000) == 0)
           disp(sprintf('Reading VPX rawdata, count %d',cnt));
       end
       %***********
       tline = fgets(fd);  % read the next line
       %***********
   end
  else
      vpx = [];
      return;
  end
  %******* when all done, throw away extra zeros in data structs
  dat = dat(1:(cnt-1),:);
  fclose(fd);
  vpx.raw = dat;
  %*****************

  %***** RUN a median filter on the data that is +/- 2 samples ********
  disp('Median filtering plus and minus 2 samples');
  mdat = dat;  % smooth traces remaining
  if (MEDFILT > 0)
    disp(sprintf('Median Filtering Eye Data (+/- %d samples)',MEDFILT));
    for k = (MEDFILT+1):(size(dat,1)-MEDFILT-1)
      mdat(k,2) = median( dat((k-MEDFILT):(k+MEDFILT),2) );
      mdat(k,3) = median( dat((k-MEDFILT):(k+MEDFILT),3) );
      mdat(k,4) = median( dat((k-MEDFILT):(k+MEDFILT),4) );
      %**********
      if (mod(k,10000) == 0)
           disp(sprintf('Med filt VPX rawdata, count %d',k));
      end
      %*********
    end
  end

  %***** upsample to 1000 hz bia linear interpolation ****************
%   disp('Up sampling 1000 hz via linear interpolation');
%   [tt2,xx2] = interp_data(mdat(:,1),mdat(:,2));
%   [tt2,yy2] = interp_data(mdat(:,1),mdat(:,3));
%   [tt2,pp2] = interp_data(mdat(:,1),mdat(:,4));
%   %***** perform Gaussian smoothing with 5ms sigma window ************
%   disp('Smoothing Gaussian Kernel (sig = 5ms)');
%   [tt3,xx3] = smooth_data(tt2,xx2,GSIG);
%   [tt3,yy3] = smooth_data(tt2,yy2,GSIG);
%   [tt3,pp3] = smooth_data(tt2,pp2,GSIG);
%   %********************************************************************
%   vpx.smo = [tt3 xx3 yy3 pp3];
  %******** for now, only median filtered for smoothed
  vpx.smo = mdat;
  disp('Completed vpx smoothing and upsampling');  

  %***************** now return to the file and build up strobe history
  [vpx.tstrobes, vpx.strobes] = read_vpx.read_vpx_strobes( filename );
  %************************************************************
  
return

function [tt2,xx2] = interp_data(tt,xx)
  %*********
  if (max(size(tt)) < 10)
      tt2 = NaN;
      xx2 = NaN;
      return;
  end
  %************
  tt2 = tt(1):0.001:(max(tt)-0.001);  % resampled data up to max time (1000 kHz)
  for zk = 1:size(tt2,2)      
      targtime = tt2(zk);
      %********** find where it is in between
      for kk = 1:(size(tt,1)-1)  % go to actual, search for where new time is located between two points
        if (targtime >= tt(kk)) & (targtime < tt(kk+1)) 
            break;
        end
      end
      %** now at kk you are in between the two points
      XA = xx(kk);   % lower value before newtime
      XB = xx(kk+1); % upper value after newtime
      t1 = targtime - tt(kk);
      t2 = tt(kk+1) - targtime;
      %*********
      val = (((t2/(t1+t2))*XA) + (t1/(t1+t2))*XB);
      xx2(zk) = val;
      %**********
  end
return;

function [tt3,xx3] = smooth_data(tt2,xx2,GSig)
      FS = 1000;   % Sampling frequency, 20khz?
      pts = GSig; % gaussian size
      dpts = (2*GSig)+1;
      posfilt = gengauss(pts,dpts); % create a gaussian filter (20 samples)
      posfilt = posfilt/sum(posfilt);
      xx3 = conv(posfilt,xx2);
      xx3 = xx3(pts:end-(pts+1)); % Window back
      xx3(1:pts) = xx3(pts+1);  % flat position on leading edge (ploblem of filtering)
      xx3((end-pts+1):end) = xx2(end-pts);  % again flat 
      tt3 = tt2;
      %*******************
return;

function x = gengauss(w,p,c)
    %GENGAUSS Generate gaussian peak of unity height.
    %   X = GENGAUSS(W)  or   X = GENGAUSS(W,P,C)
    %    where:
    %    W    is the distance from center to 2% of full height
    %    P    (optional) is the number of points x will contain
    %    C    (optional) is the point number for the peak center
    %    The peak will be centered in the vector unless P and C are specified.
    %   by Richard Kramer.
    %   Copyright (c) 1988-1993 by The MathWorks, Inc.
    %   $Revision: 1.1 $  $Date: 1993/09/01 19:37:08 $

    if nargin == 2, c = 0; end
    if nargin == 1, c = 0; p = 0; end

    if  p < 1, points = round(3 * w);
        else points = round(p);
    end

    if c < 1, center = round(points/2);
        else center = round(c);
    end

    w = w/2;
    if c < 0, center = center - c; end

    for n = 1:points,
        x(n) = exp(-((n-center)/w)^2);
    end
return;


