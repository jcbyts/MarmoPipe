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
DSI_THRESH = 0.1; % use same population as ROC analysis (DSI > 0.1)
DSI_GREEN = 0.10;
FVAL_THRESH = -10e2; % Goodness of fit
%*** setup structs to store read-in data
Abase = [];
Ubase = [];
Aamp = [];
Uamp = [];
Awid = [];
Uwid = [];
DSI = [];
%*** loop over Info files and collect desired data
xdir = dir(INFO_DATA_DIR);
N = 0;
tiny = 0.00001;
MILOCNT = 0;
 for k = 3:size(xdir,1)
       filename = [INFO_DATA_DIR,filesep,xdir(k).name];
       load(filename);
       %***** criteria for inclusion
       if (Info.DSI <= DSI_THRESH)
           continue;
       end
       
       if (Info.Tfvalstim > FVAL_THRESH)
           continue;
       end
       
       %if (Info.Afvalpre > FVAL_THRESH) || (Info.Ufvalpre > FVAL_THRESH)
       %    continue;
       %end
       
       %****** see if you can plot a normalized curve per unit
       % figure(10); hold off;
       
       
       %Info
       
       %Info.Afitpre
       %Info.Ufitpre
       
       %input('check');
       
       %if (Info.MARMO == 1)   % can exclude animal
       %    continue;
       %end
       disp(sprintf('Reading sample %d info file %s',(N+1),xdir(k).name));
       DSI = [DSI ; Info.DSI];
       %*******
       if (0)
           Abase = [Abase ; Info.Afitstim.mu(1)];
           Ubase = [Ubase ; Info.Ufitstim.mu(1)];
           Aamp = [Aamp ; Info.Afitstim.mu(2)];
           Uamp = [Uamp ; Info.Ufitstim.mu(2)];
           Awid = [Awid ; (1/(Info.Afitstim.mu(3)+tiny))];
           Uwid = [Uwid ; (1/(Info.Ufitstim.mu(3)+tiny))];
       else
           Abase = [Abase ; Info.Afitpre.mu(1)];
           Ubase = [Ubase ; Info.Ufitpre.mu(1)];
           Aamp = [Aamp ; Info.Afitpre.mu(2)];
           Uamp = [Uamp ; Info.Ufitpre.mu(2)];
           Awid = [Awid ; (1/(Info.Afitpre.mu(3)+tiny))];
           Uwid = [Uwid ; (1/(Info.Ufitpre.mu(3)+tiny))];
       end
%        if (Info.MARMO == 2)
%            MILOCNT = MILOCNT + 1;
%        end
       %*********
       N = N + 1;
       clear Info;
 end
 
%  MILOCNT

 %******* plot the results now
 H = figure;
 set(H,'Position',[100 100 1000 500]);

 %***** now perform a scatter plot *********
 %***** scatter plot stim-locked
 subplot('position',[0.1 0.2 0.225 0.45]);
 if (0)
   h = plot(Ubase,Abase,'k.'); hold on;
   set(h,'Markersize',15);
   maxo = max(max(Abase),max(Ubase)) * 1.1;
   mino = min(min(Abase),min(Ubase)) * 0.2;
   plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 else
   h = loglog(Ubase,Abase,'k.'); hold on;
   set(h,'Markersize',15);
   maxo = max(max(Abase),max(Ubase)) * 2.4;
   mino = min(min(Abase),min(Ubase)) * 0.6;
   loglog( [mino,maxo], [mino,maxo], 'k--');  % line of unity     
 end
 %labels = cellstr(num2str([1:N]'));

 axis tight;
 ax = gca;
 c = ax.FontSize;
 ax.FontSize = 11;
 V = axis;
 h = xlabel('Away');
 set(h,'FontSize',14);
 h = ylabel('Towards');
 set(h,'FontSize',14);
 set(h,'Color',[1,0,0]);
 t = title('Baseline');
 t.FontSize = 16;
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.90),sprintf('N=%d',N));
 t = text(1,150,sprintf('N=%d',N));
 t.FontSize = 12;
 %**** test for significant difference 
 pval = signrank(Ubase,Abase);
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.80),sprintf('p=%6.4f',pval));
 t = text(1,90,sprintf('p=%6.4f',pval));
 t.FontSize = 12;
 
 %************************
 subplot('position',[0.4 0.2 0.225 0.45]);
 if (0)
   h = plot(Uamp,Aamp,'k.'); hold on;
   set(h,'Markersize',15);
   %zz = find( DSI >= DSI_GREEN);
   %plot(Uamp(zz),Aamp(zz),'k.'); hold on;
   maxo = max(max(Aamp),max(Uamp)) * 1.1;
   mino = min(min(Aamp),min(Uamp)) * 0.2;
   plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 else
   h = loglog(Uamp,Aamp,'k.'); hold on;
   set(h,'Markersize',15);
   %zz = find( DSI >= DSI_GREEN);
   %plot(Uamp(zz),Aamp(zz),'k.'); hold on;
   maxo = max(max(Aamp),max(Uamp)) * 2.5;
   mino = min(min(Aamp),min(Uamp)) * 0.6;
   loglog( [mino,maxo], [mino,maxo], 'k--');  % line of unity  
 end
 %************
 axis tight;
 ax = gca;
 c = ax.FontSize;
 ax.FontSize = 11;
 V = axis;
 h = xlabel('Away');
 set(h,'FontSize',14);
 h = ylabel('Towards');
 set(h,'Color',[1,0,0]);
 set(h,'FontSize',14);
 t = title('Amplitude');
 t.FontSize = 16;
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.90),sprintf('N=%d',N));
 t = text(3,250,sprintf('N=%d',N));
 t.FontSize = 12;
 %**** test for significant difference 
 pval = signrank(Uamp,Aamp);
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.80),sprintf('p=%6.4f',pval));
 t = text(3,160,sprintf('p=%6.4f',pval));
 t.FontSize = 12;
 
 %************************
 subplot('position',[0.7 0.2 0.225 0.45]);
 h = plot(log(Uwid),log(Awid),'k.'); hold on;
 set(h,'Markersize',15);
 %zz = find( DSI >= DSI_GREEN);
 %plot(log(Uwid(zz)),log(Awid(zz)),'k.'); hold on;
 maxo = log( max(max(Awid),max(Uwid)) * 2.0);
 mino = log( min(min(Awid),min(Uwid)) * 0.2);
 plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 %**********
 axis tight;
 ax = gca;
 c = ax.FontSize;
 ax.FontSize = 11;
 V = axis;
 h = xlabel('Away');
 set(h,'FontSize',14);
 h = ylabel('Towards');
 set(h,'Color',[1,0,0]);
 set(h,'FontSize',14);
 t = title('Log Width');
 t.FontSize = 16;
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.90),sprintf('N=%d',N));
 t = text(-4,3.5,sprintf('N=%d',N));
 t.FontSize = 12;
 %**** test for significant difference 
 pval = signrank(Uwid,Awid);
 %text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.80),sprintf('p=%6.4f',pval));
 t = text(-4,2.8,sprintf('p=%6.4f',pval));
 t.FontSize = 12;
 
 