function sp = auto_channel_threshold(spdata,spinfo,FullSpread)
%******* function sp = auto_channel_threshold(spdata,spinfo,FullSpread)
%**** input:  spdata is a single data channel continuous (unfiltered) data
%****         spinfo is the info struct from open Ephys
%****         FullSpread, if set will be used as mean rate to set thresh
%****           automatically during hash sorting
%****  output:  sp - struct with thresholded spike times (one channel)

% bandpass filter between 1000 and 6000 Hz
Fs = spinfo.header.sampleRate;
[b,a] = butter(3, [(1000/Fs*2), (6000/Fs*2)]);

sp = struct('ss', [], 'clu', []);

nChannels = 1;  % size(data,1);
clustOffset = 0;
refms = 1.0;
ref_period = floor(Fs*refms/1000);  % 2 ms
threshold = -4;  % SD of data, user controlled setting
thresh_gain = 0.9; % amount to lower threshold each pass
thresh_grow = (1/0.95);  % growth in case too much rate
thresh_error = 0.20;  % get within 20% of desired
upswing = 0;  % control variable, so adjust growth for convergence (small steps)
MAXSTEPS = 10;  % maximum thresholding steps, then stop if not converged
steps = 0;
ZOOMTIME = 5;  % in seconds, for a few zoom in plots

for iCh = 1:nChannels
    
  fprintf('Sorting channel %d\n',iCh);
  
  zdata = spdata(:,iCh);
  zdata = filter(b,a,zdata);
  zdata = flipud(zdata);
  zdata = filter(b,a,zdata);
  zdata = flipud(zdata);
  SD_data = std(zdata);
  %***** for plotting pick couple zoom in times
  za(1) = floor(size(zdata,1)/4);
  zb(1) = za(1) + floor( ZOOMTIME * Fs);
  za(2) = floor(3*size(zdata,1)/4);
  zb(2) = za(2) + floor( ZOOMTIME * Fs);
  zcolo = 'bg';
  %*****************
  
  fthresh = threshold * SD_data;
  got_thresh = 0;  % cycle command line to pick threshold

  while (got_thresh == 0)

    %********* plot the total data stream
   
    if ~FullSpread
      hf = figure(1); clf
      set(hf,'position',[100 300 1200  600]);
      h1 = subplot('position',[0.05 0.6 0.9 0.25]);
      plot(zdata,'k-'); hold on;
      axis tight;
      ylabel('mV');
      xlabel('samples');
      V = axis;
      Vax = abs( fthresh * 6);  % zoom in near thresh
      axis([V(1) V(2) -Vax Vax]);
      plot([V(1),V(2)],[fthresh,fthresh],'r--');
      for ik = 1:2
        plot([za(ik),za(ik)],[-Vax,Vax],[zcolo(ik),'-']);
        plot([zb(ik),zb(ik)],[-Vax,Vax],[zcolo(ik),'-']);
      end
      %********** replot a couple zoom ins
      for ik = 1:2
        h2 = subplot('position',[0.05 0.05+(ik-1)*0.3 0.9 0.25]);
        plot(zdata(za(ik):zb(ik)),[zcolo(ik),'-']); hold on;
        axis tight;
        ylabel('mV');
        xlabel('samples');
        V = axis;
        Vax = abs( fthresh * 6);  % zoom in near thresh
        axis([V(1) V(2) -Vax Vax]);
        plot([V(1),V(2)],[fthresh,fthresh],'r--');
      end
    end
    %************* done with plotting
    
    %******** compute the spike threshold crossing and mean rate
    ss = [];
    s = find( (zdata(2:end ) < fthresh) & ...
               (zdata(1:(end - 1)) >= fthresh) );
    if (~isempty(s))
      s = s(s > 10 & s < size(zdata, 1) - 25);    % remove spikes close to boundaries
      %**** impose refractory
      z = find( ( s(2:end) - s(1:(end-1))) > ref_period);
      if (~isempty(z))
        ss = s(z);
      end
    end
    %*******
    if (~isempty(ss))
      clu = ones(size(ss));  % all same for now, just one MU unit
    else
      clu = [];
    end
    %*******************************************

    mrate = size(ss,1)/(size(zdata,1)/Fs);
    if (FullSpread)
       disp(sprintf('Mean Rate %5.2f Hz Chan(%d) thresh %5.2f',mrate,iCh,threshold));
       steps = steps + 1;
       if (mrate < FullSpread)
           threshold = threshold * thresh_gain;
           nval = threshold;
           upswing = 1;
       else
           error = (mrate - FullSpread)/FullSpread;
           if (error < thresh_error)
              nval = [];  % auto accept the thresh
           else
              if (upswing == 1) 
                thresh_gain = 0.5*(thresh_gain + 1.0);  % reduce shrinkage
                thresh_grow = (1/thresh_gain);
              end
              threshold = threshold * thresh_grow;
              nval = threshold;
              upswing = 0;
           end
       end
       if (steps > MAXSTEPS)
           nval = [];  % stop where you are, too long to converge
       end
    else    
       title(sprintf('Mean Rate %5.2f Hz Chan(%d)',mrate,iCh));
       fprintf('Current threshold: %5.2f std\n',threshold);
       nval = input('Return to accept, or enter new value: ');
    end
    %***************
    if isempty(nval)
       got_thresh = 1;
       %****** finalize spike structure
       sp.ss = [sp.ss; ss(:)];
       sp.clu = [sp.clu; clu(:)+clustOffset];
    
       clustOffset = max(sp.clu);
       %***************
    else
        threshold = nval;
        fthresh = threshold * SD_data;
    end    
  end  % end of while loop
end % end channels
disp('Finished sorting');

%% save spikes
sp.cids = unique(sp.clu)';

return;

