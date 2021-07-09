%*********** SACCADE SCRIPT *******************************
%***********************************************************
%*** NOTE:  The goal here is just to analyze a single file
%***        and do saccade flagging, plus any type of plotting
%************************************************************

if ~exist('ShowTrials', 'var')
    ShowTrials = false;  % plot trials as you go
end

%********** loop over the number of trials, find CamoFlag trials
if ShowTrials
  H = figure(1); clf
  set(H,'position',[100 200 800 800]);
end
%***** loop over trials
for i = 1:size(Exp.D,1)
    %******* check for timing strobe
    sVPX = Exp.D{i}.START_VPX;
    eVPX = Exp.D{i}.END_VPX;
    if (~isnan(sVPX) && ~isnan(eVPX) )  % make sure it was recorded OK
        %********** matlab eye position
        tt = Exp.D{i}.eyeData(:,1) - Exp.D{i}.eyeData(6,1); % sub first time
        xx = Exp.D{i}.eyeData(:,2);
        yy = Exp.D{i}.eyeData(:,3);
        %********* get VPX eye position
        zz = find( (Exp.vpx.smo(:,1) >= sVPX) & (Exp.vpx.smo(:,1) < eVPX) );
        if isempty(zz)
            [sVPX,eVPX]
            size(Exp.vpx.smo)
            disp(sprintf('Dropped trial %d, unknown cause',i));
            Exp.D{i}.eyeSmo = [];  % if not recorded, put empty record
            Exp.D{i}.slist = []; 
            continue;
        end
        vtt = Exp.vpx.smo(zz,1) - Exp.vpx.smo(zz(1),1);  %subtract first time
        vxx = Exp.vpx.smo(zz,2);
        vyy = Exp.vpx.smo(zz,3);
        vpp = Exp.vpx.smo(zz,4);
        %**** eye transform
        pixPerDeg = Exp.S.pixPerDeg;
        c = Exp.D{i}.c;
        dx = Exp.D{i}.dx;
        dy = Exp.D{i}.dy;
        %******transform positions
        xx = (xx - c(1))/(dx*pixPerDeg);
        yy = (yy - c(2))/(dy*pixPerDeg);
        %*** same on VPX data
        vxx = (vxx - c(1))/(dx*pixPerDeg);
        vyy = 1 - vyy;
        vyy = (vyy - c(2))/(dy*pixPerDeg);
        %******** process eye position to smooth and comp vel
        eye = [vtt,vxx,vyy,vpp];
        eyeSmo = +saccadeflag.eyeposition_smooth(eye);
        Exp.D{i}.eyeSmo = eyeSmo;
        %******** process smoothed eye to flag saccades
        slist = +saccadeflag.flag_saccades(eyeSmo);
        %********* store new saccade list
        Exp.D{i}.slist = slist;
        %*********************************************
        
        %********* plot eye pos over time
        if ShowTrials
            subplot('position',[0.075 0.575 0.85 0.35]); hold off;
            plot(tt,xx,'r:'); hold on;
            plot(tt,yy,'b:'); hold on;
            legend('H','V');
            plot(vtt,vxx,'r.'); hold on;
            plot(vtt,vyy,'b.'); hold on;
            %**** plot smoothed points
            plot(eyeSmo(:,1),eyeSmo(:,2),'r-');
            plot(eyeSmo(:,1),eyeSmo(:,3),'b-');
            V = axis;
            if (~isempty(slist))
               for k = 1:size(slist,1)
                   if (slist(k,7) == 1)
                       colo = 'g';
                   else
                       if (slist(k,7) == 4)
                           colo = 'm';
                       else
                           colo = 'k';
                       end
                   end
                   plot([slist(k,1),slist(k,1)],[V(3),V(4)],[colo,'-']);
                   plot([slist(k,2),slist(k,2)],[V(3),V(4)],[colo,'-']);
               end
            end
            %*******
            title(sprintf('Trial = %d',i));
            xlabel('Time (ms)');
            ylabel('EyePos');
            %*******
            if isfield(Exp.D{i}.PR,'error')
                Error = Exp.D{i}.PR.error;
            else
                Error = 0;
            end
            %****************
            subplot('position',[0.1 0.075 0.4 0.4]); hold off;
            plot(vxx,vyy,'k.'); hold on;
            plot(eyeSmo(:,2),eyeSmo(:,3),'k-');
            plot(0,0,'k+');
            axis([-10 10 -10 10]);
            title(sprintf('Error = %d',Error));
            xlabel('X Pos)');
            ylabel('Y Pos');    
        end 
        %****************
        if ShowTrials
           input('check');
        else
           fprintf('Processing trial %d\n',i);
        end
        %************
    else
        Exp.D{i}.eyeSmo = [];  % if not recorded, put empty record
        Exp.D{i}.slist = []; 
    end
end  % loop over all trials
