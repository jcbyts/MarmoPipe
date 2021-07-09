%**** POPULATION SCRIPT 
%****** runs analysis across the whole population
%****** and any stats needed

%% ***** get files from Info folder *****
athome = 0;
if (athome == 1)
   INFO_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\Info4';
else
   %INFO_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\Info';
   %INFO_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\Info2';
   INFO_DATA_DIR = 'C:\Users\Shanna\Dropbox\FovealTrans\DataAnalysis\Info4';
end

%********* Define some analysis Windows
PreWin = [-100,0];
DSI_THRESH = 0.1;
%*** setup structs to store read-in data
PopARoc = [];
PopURoc = [];
PopTag = [];
%*******
ScatARoc = [];
ScatURoc = [];
%*** loop over Info files and collect desired data
xdir = dir(INFO_DATA_DIR);
N = 0;
for k = 3:size(xdir,1)
       filename = [INFO_DATA_DIR,filesep,xdir(k).name];
       %disp(sprintf('Reading info file %s',filename));
       load(filename);
       %***** criteria for inclusion
       if (Info.DSI <= DSI_THRESH)
           continue;
       end
       %if (Info.MARMO == 2)   % can exclude animal
       %    continue;
       %end
       %*********
       tag = 0;
       %if (strcmp(Info.FileTag,'Ellie_120519'))
       %if (strcmp(Info.FileTag,'Ellie_220519'))
       %if (strcmp(Info.FileTag,'Ellie_130519'))
       % if (strcmp(Info.FileTag,'Ellie_310119'))
       %if (strcmp(Info.FileTag,'Milo_170619'))    % nice for AUC, but MU
       %if (strcmp(Info.FileTag,'Milo_190919'))     % nice but few trials
       if (strcmp(Info.FileTag,'Ellie_100219'))
         if (Info.SPClust == 1)
              tag = 1;
              Info.FileTag
              Info.SPClust
              Info.Iso
           end
       end    
       %********
       disp(sprintf('Reading info file %s',filename));
       %*******
       NPreBin = Info.NPreBin;
       TimeRoc = Info.TimeRoc;  % all the same?
       xx = 1:(NPreBin-1);
       PopARoc = [PopARoc ; Info.ARoc(xx,1)'];
       PopURoc = [PopURoc ; Info.URoc(xx,1)'];
       PopTag = [PopTag ; tag];
       PopTime = TimeRoc(xx);
       %****** compute scat stats here
       %ScatARoc = [ScatARoc ; Info.ARoc(end)];
       %ScatURoc = [ScatURoc ; Info.URoc(end)];
       %********
       StimTime = Info.TimeRoc; 
       zz = find( (StimTime(xx) >= PreWin(1)) & (StimTime(xx) < PreWin(2)) );
       ascat = mean( Info.ARoc(xx(zz),1) );
       uscat = mean( Info.URoc(xx(zz),1) );
       ascatL = mean( Info.ARoc(xx(zz),2) );
       uscatL = mean( Info.URoc(xx(zz),2) );
       ascatH = mean( Info.ARoc(xx(zz),3) );
       uscatH = mean( Info.URoc(xx(zz),3) );
       %****** transform to 1 sem
       ascatLz = ((3*ascatL) + ascatH)/4;
       ascatHz = (ascatL + (3*ascatH))/4;
       uscatLz = ((3*uscatL) + uscatH)/4;
       uscatHz = (uscatL + (3*uscatH))/4;
       %************
       ScatARoc = [ScatARoc ; [ascat,ascatLz,ascatHz]];
       ScatURoc = [ScatURoc ; [uscat,uscatLz,uscatHz]];
       
       if (tag == 1)
          figure(10);
          plot(xx,Info.ARoc(xx,1),'r-'); hold on;
          plot(xx,Info.URoc(xx,1),'b-');
          Info.FileTag
          Info.SPClust
          [ascat,uscat]
          input('stop');
          
       end
       
       %***** debug and check
       if (0)
         figure(10); hold off;
         plot(StimTime(xx),Info.ARoc(xx,1),'r:'); hold on;
         plot(StimTime(xx(zz)),Info.ARoc(xx(zz),1),'r-'); hold on;
         plot(StimTime(xx),Info.URoc(xx,1),'k:'); hold on;
         plot(StimTime(xx(zz)),Info.URoc(xx(zz),1),'k-'); hold on;
         xlabel('Time');
         ylabel('AUC');
         title(sprintf('Unit(%s) (%4.2f,%4.2f)',xdir(k).name,ascat,uscat));
         input('check');
       end
       %*********
       N = N + 1;
       clear Info;
end
%******* compute the mean stats for average normalized PSTH
AuuRoc = mean(PopARoc);
AsuRoc = std(PopARoc)/sqrt(N);
UuuRoc = mean(PopURoc);
UsuRoc = std(PopURoc)/sqrt(N);

 %******* plot the results now
 H = figure;
 set(H,'Position',[100 100 1200 600]);
 %****** general plot params
 SM = 2;

 %******* plot sac locked PSTH
 subplot('position',[0.1 0.2 0.4 0.6]);
 %****** subset of trials where no ori change occured
 %****** subset of trials where no ori change occured
 h2 = plot(PopTime,AuuRoc,'r-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(PopTime,AuuRoc+(SM*AsuRoc),'r-');hold on;
 h2b = plot(PopTime,AuuRoc-(SM*AsuRoc),'r-');hold on;
 %****** subset of trials where ori change DID occur
 h2 = plot(PopTime,UuuRoc,'b-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(PopTime,UuuRoc+(SM*UsuRoc),'b-');hold on;
 h2b = plot(PopTime,UuuRoc-(SM*UsuRoc),'b-');hold on;
 %************
 axis tight;
 V = axis;
 Vmax = 1;
 Vmin = 0.4;
 %axis([PopTime(1) PopTime(end) Vmin Vmax]);
 axis([-200 PopTime(end) Vmin Vmax]);
 ax = gca;
 c = ax.FontSize;
 ax.FontSize = 16;
 plot([0,0],[Vmin,Vmax],'k-');
 plot([PopTime(1),PopTime(end)],[0.5,0.5],'k--');
 plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
 plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
 %axis(FontSize,14)
 t = xlabel('Time (ms)');
 t(1).FontSize = 18;
 t = ylabel('AUC');
 t(1).FontSize = 18;
 t = title('Discrimation');
 t(1).FontSize = 22;
 %text((V(1)+(V(2)-V(1))*0.2),(V(3)+(V(4)-V(3))*0.7),sprintf('N=%d',N));
 t = text(-175,0.9,sprintf('N=%d',N));
 t(1).FontSize = 18;
 
 %***** now perform a scatter plot *********
 %***** scatter plot stim-locked
 subplot('position',[0.6 0.2 0.3 0.6]);
 for k = 1:size(ScatURoc,1)
     colo = 'k';
     if (PopTag(k) == 1)
         continue; %colo = 'r';
     end
     h = plot(ScatURoc(k,1),ScatARoc(k,1),[colo,'.']); hold on;
     set(h,'Markersize',20);
     plot([ScatURoc(k,1),ScatURoc(k,1)],...
          [ScatARoc(k,2),ScatARoc(k,3)],[colo,'-']); hold on;
     plot([ScatURoc(k,2),ScatURoc(k,3)],...
          [ScatARoc(k,1),ScatARoc(k,1)],[colo,'-']); hold on;
 end
 for k = 1:size(ScatURoc,1)
     colo = 'k';
     if (PopTag(k) == 0)
         continue; 
     else
         colo = 'k';
     end
     h = plot(ScatURoc(k,1),ScatARoc(k,1),[colo,'.']); hold on;
     set(h,'Markersize',20);
     plot([ScatURoc(k,1),ScatURoc(k,1)],...
          [ScatARoc(k,2),ScatARoc(k,3)],[colo,'-']); hold on;
     plot([ScatURoc(k,2),ScatURoc(k,3)],...
          [ScatARoc(k,1),ScatARoc(k,1)],[colo,'-']); hold on;
  end 
%  labels = cellstr(num2str([1:N]'));
%  text(ScatARoc,ScatURoc,labels,'VerticalAlignment','bottom',...
%                                'HorizontalAlignment','right');
 maxo = 1;
 mino = 0.4;
 plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 axis tight;
 ax = gca;
 c = ax.FontSize;
 ax.FontSize = 16;
 V = axis;
 t = xlabel('Saccade Away Sensitivity');
 t(1).FontSize = 18;
 t(1).Color = 'b';
 t = ylabel('Saccade Towards Sensitivity');
 t(1).FontSize = 18;
 t(1).Color = 'r';
 t = title('Population Sensitivity');
 t(1).FontSize = 22;
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.8),sprintf('N=%d',N));
 t = text(0.45,0.9,sprintf('N=%d',N));
 t(1).FontSize = 16;
 %**** test for significant difference 
 pval = signrank(ScatURoc(:,1),ScatARoc(:,1));
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.6),sprintf('p=%6.4f',pval));
 t = text(0.45,0.86,sprintf('p=%6.4f',pval));
 t(1).FontSize = 16;
 
 
 %*********
 H = figure;
 set(H,'Position',[200 200 1200 600]);
 %***** now perform a scatter plot *********
 %***** scatter plot stim-locked
 subplot('position',[0.6 0.2 0.3 0.6]);
 for k = 1:size(ScatURoc,1)
     colo = 'w';
     if (PopTag(k) == 1)
         colo = 'k';
     end
     h = plot(ScatURoc(k,1),ScatARoc(k,1),[colo,'.']); hold on;
     if (PopTag(k) == 1)
         set(h,'Markersize',20);
     else
         set(h,'Markersize',1);
     end
     if (PopTag(k) == 1)
       plot([ScatURoc(k,1),ScatURoc(k,1)],...
          [ScatARoc(k,2),ScatARoc(k,3)],[colo,'-']); hold on;
       plot([ScatURoc(k,2),ScatURoc(k,3)],...
          [ScatARoc(k,1),ScatARoc(k,1)],[colo,'-']); hold on;
     end
 end
%  labels = cellstr(num2str([1:N]'));
%  text(ScatARoc,ScatURoc,labels,'VerticalAlignment','bottom',...
%                                'HorizontalAlignment','right');
 maxo = 1;
 mino = 0.4;
 plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 axis tight;
 ax = gca;
 c = ax.FontSize;
 ax.FontSize = 16;
 V = axis;
 t = xlabel('Saccade Away Sensitivity');
 t(1).FontSize = 18;
 t(1).Color = 'b';
 t = ylabel('Saccade Towards Sensitivity');
 t(1).FontSize = 18;
 t(1).Color = 'r';

 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.8),sprintf('N=%d',N));
%  t = text(0.45,0.9,sprintf('N=%d',N));
%  t(1).FontSize = 16;
 %**** test for significant difference 
 pval = signrank(ScatURoc(:,1),ScatARoc(:,1));
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.6),sprintf('p=%6.4f',pval));
%  t = text(0.45,0.86,sprintf('p=%6.4f',pval));
%  t(1).FontSize = 16;