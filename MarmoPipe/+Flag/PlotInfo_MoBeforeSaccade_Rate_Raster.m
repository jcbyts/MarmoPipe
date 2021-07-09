function PlotInfo_MoBeforeSaccade_Rate_Raster(Info)
  
  %******** download variables stored in info
  fields = fieldnames(Info);
  for k = 1:size(fields,1)
      str = [fields{k} ' = Info.' fields{k} ';'];
      eval(str);
  end
  %*******************************
  
  
  maxo = max(max(AuuStim+(2*AsuStim)),max(UuuStim+(2*UsuStim)));
  maxo = maxo * 1.1;
  MR = 4;
  synchscale = 0;
  
  %****** plot the results ****************
  H = figure;
  set(H,'Position',[100 100 400 600]);
  
  
  %**************************************
  %***************************************  
  subplot('position',[0.2 0.50 0.7 0.45]);
  tick = 2;
  MR = 3;
  %********
  
  
  %**********
  if (0) 
    zb = find( ChangList == 1);
    zbb = ismember(StimRastList(:,2),zb);
    zbb2 = find( zbb == 1 );
    if (1)
     %hh = plot(StimRastList(zbb,1),StimRastList(zbb,2),'ro'); hold on;
     %set(hh,'Markersize',MR);
     
     hh = plot(StimRastList(zbb2,1),StimRastList(zbb2,2),'r.'); hold on;
     set(hh,'Markersize',MR);
    else  
     for ii = 1:size(zbb2,1)
        hh = plot([StimRastList(zbb2(ii),1),StimRastList(zbb2(ii),1)],...
                  [(StimRastList(zbb2(ii),2)-tick),(StimRastList(zbb2(ii),2)+tick)],'r-'); hold on;
        set(hh,'LineWidth',1);
     end
    end
    %******
    zc = find( ChangList == 2);  % blank trials
    zcc = ismember(StimRastList(:,2),zc);
    zcc2 = find( zcc == 1);
    if (1)
     %hh = plot(StimRastList(zcc,1),StimRastList(zcc,2),'bo'); hold on;
     %set(hh,'Markersize',MR);
     
     hh = plot(StimRastList(zcc2,1),StimRastList(zcc2,2),'b.'); hold on;
     set(hh,'Markersize',MR); 
    else
     for ii = 1:size(zcc2,1)     
        hh = plot([StimRastList(zcc2(ii),1),StimRastList(zcc2(ii),1)],...
                  [(StimRastList(zcc2(ii),2)-tick),(StimRastList(zcc2(ii),2)+tick)],'b-'); hold on;
        set(hh,'LineWidth',1);
     end
    end
  else
    imo = ones(SacNum,(BefStim+200),3);
    %******* red points
    zb = find( ChangList == 1);
    zbb = ismember(StimRastList(:,2),zb);
    zbb2 = find( zbb == 1 );
    for ii = 1:size(zbb2,1)
        if (StimRastList(zbb2(ii),1) < 200)
          imo(StimRastList(zbb2(ii),2),floor(BefStim+StimRastList(zbb2(ii),1)),2:3) = 0;
          if (StimRastList(zbb2(ii),2) < SacNum)
             imo(StimRastList(zbb2(ii),2)+1,floor(BefStim+StimRastList(zbb2(ii),1)),2:3) = 0;
          end
          if (StimRastList(zbb2(ii),2) > 1)
             imo(StimRastList(zbb2(ii),2)-1,floor(BefStim+StimRastList(zbb2(ii),1)),2:3) = 0;
          end
        end
    end
    %******* red points
    zc = find( ChangList == 2);
    zcc = ismember(StimRastList(:,2),zc);
    zcc2 = find( zcc == 1 );
    for ii = 1:size(zcc2,1)
        if (StimRastList(zcc2(ii),1) < 200)
           imo(StimRastList(zcc2(ii),2),floor(BefStim+StimRastList(zcc2(ii),1)),1:2) = 0;
           if (StimRastList(zcc2(ii),2) < SacNum)
             imo(StimRastList(zcc2(ii),2)+1,floor(BefStim+StimRastList(zcc2(ii),1)),1:2) = 0;
           end
           if (StimRastList(zcc2(ii),2) > 1)
             imo(StimRastList(zcc2(ii),2)-1,floor(BefStim+StimRastList(zcc2(ii),1)),1:2) = 0;
           end
        end
    end
    %*****************
    imagesc(-BefStim:200,0:SacNum,flipud(imo)); hold on;
  end
  %********
  plot([StimWin(1),StimWin(1)],[0,SacNum],'k--');
  plot([StimWin(2),StimWin(2)],[0,SacNum],'k--');
  %*******
  %axis([-BefStim AftStim 0 SacNum]);
  axis([-BefStim 200 0 SacNum]);
  h = plot([0,0],[0,SacNum],'k-');
  set(h,'Linewidth',2);
  %h = xlabel('Time (ms)');
  %set(h,'FontSize',14);
  h = ylabel('Trials');
  set(h,'FontSize',14);
  h = title('From Stimulus Onset'); 
  set(h,'FontSize',14);
  
  %***** plot the psth *******************
  hold on;
  subplot('position',[0.2 0.1 0.7 0.30]);
  %****** subset of trials where no ori change occured
  h2 = plot(-BefStim:AftStim,AuuStim,'r-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(-BefStim:AftStim,AuuStim+(2*AsuStim),'r-');hold on;
  h2b = plot(-BefStim:AftStim,AuuStim-(2*AsuStim),'r-');hold on;
  %****** subset of trials where ori change DID occur
  h2 = plot(-BefStim:AftStim,UuuStim,'b-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(-BefStim:AftStim,UuuStim+(2*UsuStim),'b-');hold on;
  h2b = plot(-BefStim:AftStim,UuuStim-(2*UsuStim),'b-');hold on;
  %************
  axis tight;
  V = axis;
  if (synchscale)
     Vmax = maxo;
  else
      Vmax = V(4);
  end
  
  Vmax = 40;
  %axis([-BefStim AftStim 0 Vmax]);
  axis([-BefStim 200 0 Vmax]);
  h = plot([0,0],[0,Vmax],'k-');
  set(h,'Linewidth',2);
  plot([StimWin(1),StimWin(1)],[0,Vmax],'k--');
  plot([StimWin(2),StimWin(2)],[0,Vmax],'k--');
  h = xlabel('Time (ms)');
  set(h,'FontSize',14);
  h = ylabel('Rate (sp/s)');
  set(h,'FontSize',14);
  h = title('PSTH');
  set(h,'FontSize',14);
  
  H = figure;
  set(H,'Position',[600 100 400 600]);
  
  subplot('position',[0.2 0.50 0.7 0.45]);
  %********
  if (0)
    zb = find( ChangList == 1);
    zbb = ismember(SacRastList(:,2),zb);
    hh = plot(SacRastList(zbb,1),SacRastList(zbb,2),'r.'); hold on;
    set(hh,'Markersize',MR);
    %******
    zc = find( ChangList == 2);  % blank trials
    zcc = ismember(SacRastList(:,2),zc);
    hh = plot(SacRastList(zcc,1),SacRastList(zcc,2),'b.'); hold on;
    set(hh,'Markersize',MR);
    plot([PreWin(1),PreWin(1)],[0,SacNum],'k--');
    plot([PreWin(2),PreWin(2)],[0,SacNum],'k--');
  else
    imo = ones(SacNum,(BefSac+50),3);
    %******* red points
    zb = find( ChangList == 1);
    zbb = ismember(SacRastList(:,2),zb);
    zbb2 = find( zbb == 1 );
    for ii = 1:size(zbb2,1)
        if (SacRastList(zbb2(ii),1) < 200)
          imo(SacRastList(zbb2(ii),2),floor(BefSac+SacRastList(zbb2(ii),1)),2:3) = 0;
          if (SacRastList(zbb2(ii),2) < SacNum)
             imo(SacRastList(zbb2(ii),2)+1,floor(BefSac+SacRastList(zbb2(ii),1)),2:3) = 0;
          end
          if (SacRastList(zbb2(ii),2) > 1)
             imo(SacRastList(zbb2(ii),2)-1,floor(BefSac+SacRastList(zbb2(ii),1)),2:3) = 0;
          end
%           if (SacRastList(zbb2(ii),2) < SacNum-1)
%              imo(SacRastList(zbb2(ii),2)+2,floor(BefSac+SacRastList(zbb2(ii),1)),2:3) = 0;
%           end
%           if (SacRastList(zbb2(ii),2) > 2)
%              imo(SacRastList(zbb2(ii),2)-2,floor(BefSac+SacRastList(zbb2(ii),1)),2:3) = 0;
%           end
        end
    end
    
    %******* red points
    zc = find( ChangList == 2);
    zcc = ismember(SacRastList(:,2),zc);
    zcc2 = find( zcc == 1 );
    for ii = 1:size(zcc2,1)
        if (SacRastList(zcc2(ii),1) < 200)
           imo(SacRastList(zcc2(ii),2),floor(BefSac+SacRastList(zcc2(ii),1)),1:2) = 0;
           if (SacRastList(zcc2(ii),2) < SacNum)
             imo(SacRastList(zcc2(ii),2)+1,floor(BefSac+SacRastList(zcc2(ii),1)),1:2) = 0;
           end
           if (SacRastList(zcc2(ii),2) > 1)
             imo(SacRastList(zcc2(ii),2)-1,floor(BefSac+SacRastList(zcc2(ii),1)),1:2) = 0;
           end
%            if (SacRastList(zcc2(ii),2) < SacNum-1)
%              imo(SacRastList(zcc2(ii),2)+2,floor(BefSac+SacRastList(zcc2(ii),1)),1:2) = 0;
%            end
%            if (SacRastList(zcc2(ii),2) > 2)
%              imo(SacRastList(zcc2(ii),2)-2,floor(BefSac+SacRastList(zcc2(ii),1)),1:2) = 0;
%            end
        end
    end
    %*****************
    imagesc(-BefSac:50,0:SacNum,flipud(imo)); hold on;      
  end
  %***********
  %axis([-BefSac AftSac 0 SacNum]);
  axis([-200 50 0 SacNum]);
  h = plot([0,0],[0,SacNum],'k-');
  set(h,'Linewidth',2);
  %h = xlabel('Time (ms)');
  %set(h,'FontSize',14);
  h = ylabel('Trials');
  set(h,'FontSize',14);
  h = title('From Saccade Onset');   
  set(h,'FontSize',14);
  
  %***** plot the psth *******************
  hold on;
  subplot('position',[0.2 0.1 0.7 0.30])
  %****** subset of trials where no ori change occured
  h2 = plot(-BefSac:AftSac,AuuPre,'r-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(-BefSac:AftSac,AuuPre+(2*AsuPre),'r-');hold on;
  h2b = plot(-BefSac:AftSac,AuuPre-(2*AsuPre),'r-');hold on;
  %****** subset of trials where ori change DID occur
  h2 = plot(-BefSac:AftSac,UuuPre,'b-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(-BefSac:AftSac,UuuPre+(2*UsuPre),'b-');hold on;
  h2b = plot(-BefSac:AftSac,UuuPre-(2*UsuPre),'b-');hold on;
  %************
  axis tight;
  V = axis;
  if (synchscale)
     Vmax = maxo;
  else
      Vmax = V(4);
  end
  Vmax = 40;
  %axis([-BefSac AftSac 0 Vmax]);
  axis([-200 50 0 Vmax]);
  h = plot([0,0],[0,Vmax],'k-');
  set(h,'Linewidth',2);
  plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
  plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
  h = xlabel('Time (ms)');
  set(h,'FontSize',14);
  h = ylabel('Rate (sp/s)');
  set(h,'FontSize',14);
  h = title('PSTH');
  set(h,'FontSize',14);
  
   
  return;
  
  
 function PlotWithVonMises(OriVals,TrA,zPreOri,StimSpk,Atune,colo,XJit,YJit)
  
   %******
   OriNum = size(OriVals,1);
   for i = 1:OriNum
        %****** plot raw counts
        % first plot cases where saccade to RF
        pzz = find( (zPreOri(TrA) == OriVals(i)));
        xjit = (rand(size(pzz))-0.5)*XJit;
        yjit = (rand(size(pzz))-0.5)*YJit;     
        %********
        H = plot(OriVals(i)*ones(size(pzz))+xjit,StimSpk(TrA(pzz))+yjit,[colo,'.']); hold on;
        set(H,'Markersize',10);
        %****** plot vonMises fit
        h2 = plot(OriVals,Atune.mu,[colo,'-']);
        set(h2,'Linewidth',2);
        h2 = plot(OriVals,Atune.mu + (2 * Atune.sem),[colo,'-']);
        h2 = plot(OriVals,Atune.mu - (2 * Atune.sem),[colo,'-']);
        %*********
  end
  axis tight;
  V = axis;
  axis([0 360 0 V(4)]);
  ax = gca;
  set(ax,'XTickLabel',[]);
  
return;