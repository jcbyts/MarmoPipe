function PlotInfo_MoBeforeSaccade(Info)
  
  %******** download variables stored in info
  fields = fieldnames(Info);
  for k = 1:size(fields,1)
      str = [fields{k} ' = Info.' fields{k} ';'];
      eval(str);
  end
  %*******************************
  
  %****** plot the results ****************
  H = figure;
  %set(H,'Position',[100 100 1200 800]);
  
% %   %*********** plot the RF
%   KN = size(RFsets,2);
%   for it = 1:(KN-1)
%     subplot('position',[(0.05 + 0.1*(it-1)),0.85,0.08,0.12]);
%     svec = RFsets{1,it};
%     imagesc(Zx,Zy,svec,[RFmino RFmaxo]); hold on;
%     plot([-11,11],[0,0],'k-');
%     plot([0,0],[-11,11],'k-');
%     axis off;
%     h = title(sprintf('%4.1f',-tXX(it)));
%   end
%   subplot('position',[(0.05 + 0.1*(KN-1)),0.85,0.08,0.12]);
%   imagesc(Zx,Zy,ones(size(svec))*RFmino,[RFmino RFmaxo]); hold on;
%   axis off;
%   h = colorbar;
%   
  %******** plot the motion tuning with DSI index
%   subplot('position',[0.05 0.625 0.10 0.15]);
%   px = []; py = [];
%   rx = []; ry = [];
%   for k = 1:size(OriVals,1)
%      ori = OriVals(k)*(pi/180);
%      vc = Ttunestim.mu(k) * complex(cos(ori),sin(ori));
%      px = [px ; real(vc)];
%      py = [py ; imag(vc)];
%      vc = Trawstim(k) * complex(cos(ori),sin(ori));
%      rx = [rx ; real(vc)];
%      ry = [ry ; imag(vc)];
%   end
%   px = [px ; px(1)];
%   py = [py ; py(1)];
%   rx = [rx ; rx(1)];
%   ry = [ry ; ry(1)];
%   %**** draw the circle
%   maxo = max(Ttunestim.mu);
%   rmaxo = max(max(abs(rx)),max(abs(ry))) * 1.1;
%   plot(rx,ry,'b.'); hold on;
%   plot(px,py,'b-'); hold on;
%   ori = OriVals(PrefDir)*(pi/180);
%   Pvec = complex(cos(ori),sin(ori));
%   plot([0,real(Pvec * maxo)],[0,imag(Pvec * maxo)],'g-');
%   plot(0,0,'k+');
%   h = plot([0,real(Wvec * maxo)],[0,imag(Wvec * maxo)],'k-');
%   set(h,'Linewidth',2);
%   axis([-rmaxo rmaxo -rmaxo rmaxo]);
%   title(sprintf('DSI = %4.2f',DSI));
%%*******************************************
 %******************************************* 
 %*******************************************
  
%   %********* PreSac - Plot AUC for Att and UnAtt
%   subplot('position',[0.75 0.60 0.225 0.275])
%   %****** subset of trials where no ori change occured
%   xx = 1:(NPreBin-1);
%   h2 = plot(TimeRoc(xx),ARoc(xx,1),'r-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(TimeRoc(xx),ARoc(xx,2),'r-');hold on;
%   h2b = plot(TimeRoc(xx),ARoc(xx,3),'r-');hold on;
%   %****** subset of trials where ori change DID occur
%   h2 = plot(TimeRoc(xx),URoc(xx,1),'k-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(TimeRoc(xx),URoc(xx,2),'k-');hold on;
%   h2b = plot(TimeRoc(xx),URoc(xx,3),'k-');hold on;
%   %************
%   axis tight;
%   V = axis;
%   axis([TimeRoc(1) TimeRoc(xx(end)) V(3) V(4)]);
%   plot([0,0],[V(3),V(4)],'k-');
%   plot([TimeRoc(1),TimeRoc(xx(end))],[0.5,0.5],'k--');
%   plot([PreWin(1),PreWin(1)],[V(3),V(4)],'k--');
%   plot([PreWin(2),PreWin(2)],[V(3),V(4)],'k--');
%   xlabel('Time (ms)');
%   ylabel('AUC');
%   title('PreSac');
% %   
% 
%   
%   %*****************************************************
%   %*****************************************************
%   %*****************************************************
%   %******** find max range *************
%   Vmax = max([max(AuuStimPref+(2*AsuStimPref)),max(AuuPrePref+(2*AsuPrePref)),...
%               max(UuuStimPref+(2*UsuStimPref)),max(UuuPrePref+(2*UsuPrePref))]);
%           
%   %********* Stim - attended, pref vs non-pref
%   subplot('position',[0.2 0.6 0.225 0.275]);
%   %****** subset of trials where no ori change occured
%   h2 = plot(-BefStim:AftStim,AuuStimPref,'g-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(-BefStim:AftStim,AuuStimPref+(2*AsuStimPref),'g-');hold on;
%   h2b = plot(-BefStim:AftStim,AuuStimPref-(2*AsuStimPref),'g-');hold on;
%   %****** subset of trials where ori change DID occur
%   h2 = plot(-BefStim:AftStim,AuuStimNPref,'b-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(-BefStim:AftStim,AuuStimNPref+(2*AsuStimNPref),'b-');hold on;
%   h2b = plot(-BefStim:AftStim,AuuStimNPref-(2*AsuStimNPref),'b-');hold on;
%   %************
%   axis tight;
%   V = axis;
%   axis([-BefStim AftStim 0 Vmax]);
%   plot([0,0],[0,Vmax],'k-');
%   plot([StimWin(1),StimWin(1)],[0,Vmax],'k--');
%   plot([StimWin(2),StimWin(2)],[0,Vmax],'k--');
%   xlabel('Time (ms)');
%   ylabel('Rate (hz)');
%   title('Attended Stim');
%   
%   %********* PreSac - attended, pref vs non-pref
%   subplot('position',[0.5 0.6 0.225 0.275])
%   %****** subset of trials where no ori change occured
%   h2 = plot(-BefSac:AftSac,AuuPrePref,'g-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(-BefSac:AftSac,AuuPrePref+(2*AsuPrePref),'g-');hold on;
%   h2b = plot(-BefSac:AftSac,AuuPrePref-(2*AsuPrePref),'g-');hold on;
%   %****** subset of trials where ori change DID occur
%   h2 = plot(-BefSac:AftSac,AuuPreNPref,'b-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(-BefSac:AftSac,AuuPreNPref+(2*AsuPreNPref),'b-');hold on;
%   h2b = plot(-BefSac:AftSac,AuuPreNPref-(2*AsuPreNPref),'b-');hold on;
%   %************
%   axis tight;
%   V = axis;
%   axis([-BefSac AftSac 0 Vmax]);
%   plot([0,0],[0,Vmax],'k-');
%   plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
%   plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
%   xlabel('Time (ms)');
%   ylabel('Rate (hz)');
%   title('Attended PreSac');
%   
%   %********* Stim - attended, pref vs non-pref
%   subplot('position',[0.2 0.15 0.225 0.275])
%   %****** subset of trials where no ori change occured
%   h2 = plot(-BefStim:AftStim,UuuStimPref,'g-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(-BefStim:AftStim,UuuStimPref+(2*UsuStimPref),'g-');hold on;
%   h2b = plot(-BefStim:AftStim,UuuStimPref-(2*UsuStimPref),'g-');hold on;
%   %****** subset of trials where ori change DID occur
%   h2 = plot(-BefStim:AftStim,UuuStimNPref,'b-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(-BefStim:AftStim,UuuStimNPref+(2*UsuStimNPref),'b-');hold on;
%   h2b = plot(-BefStim:AftStim,UuuStimNPref-(2*UsuStimNPref),'b-');hold on;
%   %************
%   axis tight;
%   V = axis;
%   axis([-BefStim AftStim 0 Vmax]);
%   plot([0,0],[0,Vmax],'k-');
%   plot([StimWin(1),StimWin(1)],[0,Vmax],'k--');
%   plot([StimWin(2),StimWin(2)],[0,Vmax],'k--');
%   xlabel('Time (ms)');
%   ylabel('Rate (hz)');
%   title('Unattended Stim');
%   
%   %********* PreSac - attended, pref vs non-pref
%   subplot('position',[0.5 0.15 0.225 0.275])
%   %****** subset of trials where no ori change occured
%   h2 = plot(-BefSac:AftSac,UuuPrePref,'g-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(-BefSac:AftSac,UuuPrePref+(2*UsuPrePref),'g-');hold on;
%   h2b = plot(-BefSac:AftSac,UuuPrePref-(2*UsuPrePref),'g-');hold on;
%   %****** subset of trials where ori change DID occur
%   h2 = plot(-BefSac:AftSac,UuuPreNPref,'b-');hold on;
%   set(h2,'Linewidth',2);
%   h2b = plot(-BefSac:AftSac,UuuPreNPref+(2*UsuPreNPref),'b-');hold on;
%   h2b = plot(-BefSac:AftSac,UuuPreNPref-(2*UsuPreNPref),'b-');hold on;
%   %************
%   axis tight;
%   V = axis;
%   axis([-BefSac AftSac 0 Vmax]);
%   plot([0,0],[0,Vmax],'k-');
%   plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
%   plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
%   xlabel('Time (ms)');
%   ylabel('Rate (hz)');
%   title('Unattended PreSac');
    %*****************************************************
  %*****************************************************
  %*****************************************************
  
  
%   showgreen = 0;
%   showblack = 0;
%   %******************
%   
%   subplot('position',[0.40 0.425 0.125 0.35]);
%   %********
%   zb = find( ChangList == 1);
%   zbb = ismember(StimRastList(:,2),zb);
%   hh = plot(StimRastList(zbb,1),StimRastList(zbb,2),'r.'); hold on;
%   set(hh,'Markersize',0.2);
%   %******
%   zc = find( ChangList == 2);  % blank trials
%   zcc = ismember(StimRastList(:,2),zc);
%   hh = plot(StimRastList(zcc,1),StimRastList(zcc,2),'b.'); hold on;
%   set(hh,'Markersize',0.2);
%   %********
%   plot([StimWin(1),StimWin(1)],[0,SacNum],'k--');
%   plot([StimWin(2),StimWin(2)],[0,SacNum],'k--');
%   %******
%   if (showblack)
%      hh = plot(StimSacDur(:,1),StimSacDur(:,2),'k.'); 
%   end
%   if (showgreen)
%      hh2 = plot(StimMSacDur(:,1),StimMSacDur(:,2),'g.');
%   end
%   %*******
%   axis([-BefStim AftStim 0 SacNum]);
%   plot([0,0],[0,SacNum],'k-');
%   xlabel('Time (ms)');
%   ylabel('Trials');
%   title('(Timelock Stim Onset)'); 
  
%   %***** plot the psth *******************
%   hold on;
%   subplot('position',[0.40 0.075 0.125 0.25])
%   %****** subset of trials where no ori change occured
%   h2 = plot(-BefStim:AftStim,AuuStim,'r-');hold on;
%   h2b = plot(-BefStim:AftStim,AuuStim+(2*AsuStim),'r:');hold on;
%   h2b = plot(-BefStim:AftStim,AuuStim-(2*AsuStim),'r:');hold on;
%   %****** subset of trials where ori change DID occur
%   h2 = plot(-BefStim:AftStim,UuuStim,'b-');hold on;
%   h2b = plot(-BefStim:AftStim,UuuStim+(2*UsuStim),'b:');hold on;
%   h2b = plot(-BefStim:AftStim,UuuStim-(2*UsuStim),'b:');hold on;
%   %************
%   axis tight;
%   V = axis;
%   axis([-BefStim AftStim 0 V(4)]);
%   plot([0,0],[0,V(4)],'k-');
%   plot([StimWin(1),StimWin(1)],[0,V(4)],'k--');
%   plot([StimWin(2),StimWin(2)],[0,V(4)],'k--');
%   xlabel('Time (ms)');
%   ylabel('Rate (hz)');
%   title('PSTH');
%   
%   
%   subplot('position',[0.575 0.425 0.125 0.35]);
%   %********
%   zb = find( ChangList == 1);
%   zbb = ismember(SacRastList(:,2),zb);
%   hh = plot(SacRastList(zbb,1),SacRastList(zbb,2),'r.'); hold on;
%   set(hh,'Markersize',0.2);
%   %******
%   zc = find( ChangList == 2);  % blank trials
%   zcc = ismember(SacRastList(:,2),zc);
%   hh = plot(SacRastList(zcc,1),SacRastList(zcc,2),'b.'); hold on;
%   set(hh,'Markersize',0.2);
%   plot([PreWin(1),PreWin(1)],[0,SacNum],'k--');
%   plot([PreWin(2),PreWin(2)],[0,SacNum],'k--');
%   %******
%   if (showblack)
%      hh = plot(SacSacDur(:,1),SacSacDur(:,2),'k.');      % leave in vpx
%   end
%   if (showgreen)
%      hh2 = plot(SacMSacDur(:,1),SacMSacDur(:,2),'g.');   % leave fix in matlab
%   end
%   %***********
%   axis([-BefSac AftSac 0 SacNum]);
%   plot([0,0],[0,SacNum],'k-');
%   xlabel('Time (ms)');
%   ylabel('Trials');
%   title('(Timelock Onset)');   
%   
%   %***** plot the psth *******************
%   hold on;
%   subplot('position',[0.575 0.075 0.125 0.25])
%   %****** subset of trials where no ori change occured
%   h2 = plot(-BefSac:AftSac,AuuPre,'r-');hold on;
%   h2b = plot(-BefSac:AftSac,AuuPre+(2*AsuPre),'r:');hold on;
%   h2b = plot(-BefSac:AftSac,AuuPre-(2*AsuPre),'r:');hold on;
%   %****** subset of trials where ori change DID occur
%   h2 = plot(-BefSac:AftSac,UuuPre,'b-');hold on;
%   h2b = plot(-BefSac:AftSac,UuuPre+(2*UsuPre),'b:');hold on;
%   h2b = plot(-BefSac:AftSac,UuuPre-(2*UsuPre),'b:');hold on;
%   %************
%   axis tight;
%   V = axis;
%   axis([-BefSac AftSac 0 V(4)]);
%   plot([0,0],[0,V(4)],'k-');
%   plot([PreWin(1),PreWin(1)],[0,V(4)],'k--');
%   plot([PreWin(2),PreWin(2)],[0,V(4)],'k--');
%   xlabel('Time (ms)');
%   ylabel('Rate (hz)');
%   title('PSTH');
%   
%   %******** we need some plot of saccade end points
  showrawtrace = 1;
%  subplot('position',[0.75 0.70 0.19 0.28]);
  colo = 'rkbmgc';
  for k = 1:size(TPosList,1)
     tg = TPosList(k,1);
     if (showrawtrace)
        h = plot(TPosList(k,2),TPosList(k,3),[colo(tg),'.']); hold on;
       set(h,'Markersize',4);
       %*** plot eye position before the saccade
        h = plot(TPosTrace{k}(:,2),TPosTrace{k}(:,3),[colo(tg),':']); hold on;
     end
     %******
  end
  %****** get mean eye position before saccade
  for k = 1:size(TPosTraceList,2)
      h = plot(TuPos(k,1),TuPos(k,2),[colo(k),'o']); 
      set(h,'Markersize',10); hold on;
      h = plot(TuPos(k,3),TuPos(k,4),[colo(k),'s']); 
      set(h,'Markersize',10);
  end
  %**************
  h = plot(0,0,'k+');
  set(h,'Markersize',10);
  maxo = max(max(abs(TPosList(:,2:3))));
  axis([-maxo maxo -maxo maxo]);
  ax = gca;
  c = ax.FontSize;
  ax.FontSize = 12;
%    
% Zoomed in Version
  H = figure;
  showrawtrace = 1;
%  subplot('position',[0.75 0.70 0.19 0.28]);
  colo = 'rkbmgc';
  for k = 1:size(TPosList,1)
     tg = TPosList(k,1);
     if (showrawtrace)
        h = plot(TPosList(k,2),TPosList(k,3),[colo(tg),'.']); hold on;
       set(h,'Markersize',4);
       %*** plot eye position before the saccade
        h = plot(TPosTrace{k}(:,2),TPosTrace{k}(:,3),[colo(tg),':']); hold on;
     end
     %******
  end
  %****** get mean eye position before saccade
  for k = 1:size(TPosTraceList,2)
      h = plot(TuPos(k,1),TuPos(k,2),[colo(k),'o']); 
      set(h,'Markersize',25); hold on;
      h = plot(TuPos(k,3),TuPos(k,4),[colo(k),'s']); 
      set(h,'Markersize',20);
  end
  %**************
  h = plot(0,0,'k+');
  set(h,'Markersize',10);
  maxo = max(max(abs(TPosList(:,2:3))));
  %axis([-maxo maxo -maxo maxo]);
  axis([-1 1 -1 1]);
  ax = gca;
  c = ax.FontSize;
  ax.FontSize = 12;
  set(gca,'xtick',[])
  set(gca,'ytick',[])

%*****************************************
%*****************************************
%*****************************************
%   %******** Plot stim locked tuning curve and do Von Mises Fit
%   subplot('position',[0.75 0.375 0.20 0.25]);
%   %*** plot fits over raw data points
%   zYJit = max(max(Atunestim.mu),max(Utunestim.mu)) * YJit;
%   PlotWithVonMises(OriVals,TrA,zPreOri,StimSpk{NStimBin},Atunestim,'r',XJit,zYJit);
%   PlotWithVonMises(OriVals,TrU,zPreOri,StimSpk{NStimBin},Utunestim,'b',XJit,zYJit);
%   ylabel('Spike Counts');
%   title('Stim Locked Activity');
% 
%   %*****************************************
%   subplot('position',[0.75 0.075 0.20 0.25]);
%   %*** plot fits over raw data points
%   zYJit = max(max(Atunepre.mu),max(Utunepre.mu)) * YJit;
%   PlotWithVonMises(OriVals,TrA,zPreOri,PreSpk{NPreBin},Atunepre,'r',XJit,zYJit);
%   PlotWithVonMises(OriVals,TrU,zPreOri,PreSpk{NPreBin},Utunepre,'b',XJit,zYJit);
%   ylabel('Spike Counts');
%   xlabel ('Direction')
%   title('PreSac Locked Activity');
%    
%   return;
%   
%   
%  function PlotWithVonMises(OriVals,TrA,zPreOri,StimSpk,Atune,colo,XJit,YJit)
%   
%    %******
%    OriNum = size(OriVals,1);
%    for i = 1:OriNum
%         %****** plot raw counts
%         % first plot cases where saccade to RF
%         pzz = find( (zPreOri(TrA) == OriVals(i)));
%         xjit = (rand(size(pzz))-0.5)*XJit;
%         yjit = (rand(size(pzz))-0.5)*YJit;     
%         %********
%         H = plot(OriVals(i)*ones(size(pzz))+xjit,StimSpk(TrA(pzz))+yjit,[colo,'.']); hold on;
%         set(H,'Markersize',3);
%         %****** plot vonMises fit
%         h2 = plot(OriVals,Atune.mu,[colo,'-']);
%         set(h2,'Linewidth',2);
%         h2 = plot(OriVals,Atune.mu + (2 * Atune.sem),[colo,'-']);
%         h2 = plot(OriVals,Atune.mu - (2 * Atune.sem),[colo,'-']);
%         %*********
%   end
%   axis tight;
%   V = axis;
%   axis([0 360 0 V(4)]);
%   ax = gca;
%   set(ax,'XTickLabel',[]);
  
return;