function vpx  = load_edf_file( filename )
% function  vpx = load_edf_file( filename )
% - reads edf file struct that was saved out into Matlab format
% - also detects trial start and end strobe times with taglets
%     from the fevents fields of the edf matlab struct
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

  MEDFILT = 3;   % median filtering of data series, 0 if none
  GSIG = 5;      % Gaussian smoothing in ms (must be odd integer), 0 if none
  %*******************
  disp(sprintf('Reading EDF file %s',filename));
  vpx = struct;
  vpx_strobes = struct;
  
  %****** first load the EDF file
  load(filename);
  if ~exist('edf', 'var')
      edf = load(filename);
  end
  
  if ~isstruct(edf)
     fprintf('Failed to load EDF file %s\n',filename);
     disp('Returning very unsatisfied ....');
     return;
  end
  %***** read in the raw eye data first
  MaxCnt = size(edf.FSAMPLE.px,2);
  p1 = range(edf.FSAMPLE.px(1,:));
  p2 = range(edf.FSAMPLE.px(2,:));
  EYEUSED = 2;  % how do I find this from struct ... use message?
  if (p1 > p2)  %other eye than default was used
      EYEUSED = 1;
  end           % option of binocular not considered here, if so would average both
  %********
  scale = (1/32768);   %apply transformation we used
  %**********
  tt = double(edf.FSAMPLE.time')/1000;
  xx = -scale * double( edf.FSAMPLE.px(EYEUSED,:))';
  yy = scale * double( edf.FSAMPLE.py(EYEUSED,:))';
  pp = double( edf.FSAMPLE.pa(EYEUSED,:))';
  vpx.raw = [tt xx yy pp];
  dat = vpx.raw;
   
  %******* COULD BE TOO SLOW GIVEN EDF samples 1000hz instead of other           
  %***** RUN a median filter on the data that is +/- 2 samples ********
  mdat = dat;  % smooth traces remaining
  if (MEDFILT > 0)
    disp(sprintf('Median Filtering Eye Data (+/- %d samples)',MEDFILT));
    for k = (MEDFILT+1):(size(dat,1)-MEDFILT-1)
      mdat(k,2) = median( dat((k-MEDFILT):(k+MEDFILT),2) );
      mdat(k,3) = median( dat((k-MEDFILT):(k+MEDFILT),3) );
      mdat(k,4) = median( dat((k-MEDFILT):(k+MEDFILT),4) );
      %**********
      if (mod(k,50000) == 0)
           disp(sprintf('Med filt EDF rawdata, count %d',k));
      end
      %*********
    end
  end
  %***** perform Gaussian smoothing with 5ms sigma window ************
  disp(sprintf('Smoothing Gaussian Kernel (sig = %5.1f ms)',GSIG));
  [tt3,xx3] = smooth_data(mdat(:,1),mdat(:,2),GSIG);
  [tt3,yy3] = smooth_data(mdat(:,1),mdat(:,3),GSIG);
  [tt3,pp3] = smooth_data(mdat(:,1),mdat(:,4),GSIG);
  %********************************************************************
  vpx.smo = [tt3 xx3 yy3 pp3];
  %******** for now, only median filtered for smoothed
  vpx.smo = mdat;
  vpx.dat = mdat;
  disp('Completed vpx smoothing and upsampling');  
  %***************** now return to the file and build up strobe history
  [vpx.tstrobes, vpx.strobes] = read_edf.read_edf_strobes( edf.FEVENT );
  %************************************************************
  clear edf;
  
return

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


