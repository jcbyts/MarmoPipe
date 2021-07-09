%**** POPULATION SCRIPT 
%****** runs analysis across the whole population
%****** and any stats needed

%% ***** get files from Info folder *****
athome = 0;
if (athome == 1)
   % INFO_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\Info';
   INFO_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\InfoSU_JM';
else
   %INFO_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\Info';
   %INFO_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\Info2';
   INFO_DATA_DIR = 'C:\Users\Shanna\Dropbox\FovealTrans\DataAnalysis\Info4';
end

%********* Define some analysis Windows
StimWin = [50,150];
PreWin = [-100,0];
MyBefStim = -50;
MyAftStim = 200;
MyBefPre = -200;
MyAftPre = 50;
%*** setup structs to store read-in data
PopAuuStim = [];
PopUuuStim = [];
PopAuuPre = [];
PopUuuPre = [];
%*******
ScatAuuStim = [];
ScatUuuStim = [];
ScatAuuPre = [];
ScatUuuPre = [];
%**********
NPopAuuStim = [];
NPopUuuStim = [];
NPopAuuPre = [];
NPopUuuPre = [];
%************
ISO = [];
MARMO = [];
%**************
NPopIso = [];
%******************
%*** loop over Info files and collect desired data
xdir = dir(INFO_DATA_DIR);
N = 0;
for k = 3:size(xdir,1)
       filename = [INFO_DATA_DIR,filesep,xdir(k).name];
       disp(sprintf('Reading info file %s',filename));
       load(filename); 
       %*********
%        ISO = [ISO ; Info.ISO];
%        MARMO = [MARMO ; Info.MARMO];
       %*******
       PopAuuStim = [PopAuuStim ; Info.AuuStim];
       PopUuuStim = [PopUuuStim ; Info.UuuStim];
       PopAuuPre =  [PopAuuPre  ; Info.AuuPre];
       PopUuuPre =  [PopUuuPre  ; Info.UuuPre];
       %*******
       MaxStim = max( [Info.AuuStim Info.UuuStim]);
       MaxPre = max( [Info.AuuPre Info.UuuPre]);
       NPopAuuStim = [NPopAuuStim ; (Info.AuuStim * (1/MaxStim))];
       NPopUuuStim = [NPopUuuStim ; (Info.UuuStim * (1/MaxStim))];
       NPopAuuPre =  [NPopAuuPre  ; (Info.AuuPre * (1/MaxPre))];
       NPopUuuPre =  [NPopUuuPre  ; (Info.UuuPre * (1/MaxPre))];
       %*******
       BefStim = Info.BefStim;
       AftStim = Info.AftStim;
       BefPre = Info.BefSac;
       AftPre = Info.AftSac;
       %****** compute scat stats here
       StimTime = -BefStim:AftStim;
       zz = find( (StimTime >= StimWin(1)) & (StimTime < StimWin(2)) );
       ScatAuuStim = [ScatAuuStim ; mean( Info.AuuStim(zz) )];
       ScatUuuStim = [ScatUuuStim ; mean( Info.UuuStim(zz) )];
       %********
       PreTime = -BefPre:AftPre;
       zz = find( (PreTime >= PreWin(1)) & (PreTime < PreWin(2)) );
       ScatAuuPre = [ScatAuuPre ; mean( Info.AuuPre(zz) )];
       ScatUuuPre = [ScatUuuPre ; mean( Info.UuuPre(zz) )];
       %*********
       N = N + 1;
       clear Info;
end

%******* compute the mean stats for average normalized PSTH
AuuStim = mean(NPopAuuStim);
AsuStim = std(NPopAuuStim)/sqrt(N);
UuuStim = mean(NPopUuuStim);
UsuStim = std(NPopUuuStim)/sqrt(N);
AuuPre = mean(NPopAuuPre);
AsuPre = std(NPopAuuPre)/sqrt(N);
UuuPre = mean(NPopUuuPre);
UsuPre = std(NPopUuuPre)/sqrt(N);

 %******* plot the results now
 H = figure;
 set(H,'Position',[100 100 600 300]);
 %****** general plot params
 SM = 2;
 
 %******* plot stim locked PSTH
 subplot('position',[0.1 0.2 0.35 0.6]);
 %*********
 xt = -BefStim:AftStim;
 zt = find( (xt >= MyBefStim) & (xt <= MyAftStim) );
 %****** subset of trials where no ori change occured
 h2 = plot(xt(zt),AuuStim(zt),'r-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(xt(zt),AuuStim(zt)+(SM*AsuStim(zt)),'r-');hold on;
 h2b = plot(xt(zt),AuuStim(zt)-(SM*AsuStim(zt)),'r-');hold on;
 %****** subset of trials where ori change DID occur
 h2 = plot(xt(zt),UuuStim(zt),'b-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(xt(zt),UuuStim(zt)+(SM*UsuStim(zt)),'b-');hold on;
 h2b = plot(xt(zt),UuuStim(zt)-(SM*UsuStim(zt)),'b-');hold on;
 %************
 axis tight;
 V = axis;
 Vmax = V(4);
 axis([MyBefStim MyAftStim 0 Vmax]);
 plot([0,0],[0,Vmax],'k-');
 plot([StimWin(1),StimWin(1)],[0,Vmax],'k--');
 plot([StimWin(2),StimWin(2)],[0,Vmax],'k--');
 h = xlabel('Time (ms)');
 set(h,'FontSize',14);
 h = ylabel('Norm. Rate');
 set(h,'FontSize',14);
 h= title('Stim-Locked Towards vs Away');
 set(h,'FontSize',14);
 h = text((V(1)+(V(2)-V(1))*0.2),(V(3)+(V(4)-V(3))*0.7),sprintf('N=%d',N));
 set(h,'FontSize',14);
 %text((V(2)*0.2),(V(4)*0.7),sprintf('N=%d',N));
 %******************
 
 %******* plot stim locked PSTH
 subplot('position',[0.55 0.2 0.35 0.6]);
 xt = -BefPre:AftPre;
 zt = find( (xt >= MyBefPre) & (xt <= MyAftPre) );
 %****** subset of trials where no ori change occured
 h2 = plot(xt(zt),AuuPre(zt),'r-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(xt(zt),AuuPre(zt)+(SM*AsuPre(zt)),'r-');hold on;
 h2b = plot(xt(zt),AuuPre(zt)-(SM*AsuPre(zt)),'r-');hold on;
 %****** subset of trials where ori change DID occur
 h2 = plot(xt(zt),UuuPre(zt),'b-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(xt(zt),UuuPre(zt)+(SM*UsuPre(zt)),'b-');hold on;
 h2b = plot(xt(zt),UuuPre(zt)-(SM*UsuPre(zt)),'b-');hold on;
 %************
 axis tight;
 V = axis;
 Vmax = V(4);
 axis([MyBefPre MyAftPre 0 Vmax]);
 plot([0,0],[0,Vmax],'k-');
 plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
 plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
 h = xlabel('Time (ms)');
 set(h,'FontSize',14);
 h = ylabel('Norm. Rate');
 set(h,'FontSize',14);
 h = title('Sac-Locked Towards vs Away');
 %text((V(2)*0.2),(V(4)*0.7),sprintf('N=%d',N));
 h = text((V(1)+(V(2)-V(1))*0.2),(V(3)+(V(4)-V(3))*0.7),sprintf('N=%d',N));
 set(h,'FontSize',14);
 
 %***** subsets
 amu = find( (ISO < 3) & (MARMO == 1) );
 aio = find( (ISO >= 3) & (MARMO == 1) );
 bmu = find( (ISO < 3) & (MARMO == 2) );
 bio = find( (ISO >= 3) & (MARMO == 2) );
 %**************
 
 %***** now perform a scatter plot *********
 %******* plot the results now
 H = figure;
 set(H,'Position',[100 100 1200 600]);
 %***** scatter plot stim-locked
 subplot('position',[0.2 0.2 0.3 0.6]);
 plot(ScatAuuStim,ScatUuuStim,'k.'); 
 maxo = max(max(ScatAuuStim),max(ScatUuuStim)) * 1.2;
 mino = min(min(ScatAuuStim),min(ScatUuuStim)) * 0.8;
 if (1)  % show all the same
   h = loglog( ScatUuuStim, ScatAuuStim, 'k.'); hold on;
   set(h,'Markersize',18);
 else
   h = loglog( ScatUuuStim(amu), ScatAuuStim(amu), 'ko'); hold on;
   set(h,'Markersize',5);
   h = loglog( ScatUuuStim(aio), ScatAuuStim(aio), 'k.'); hold on;
   set(h,'Markersize',15);
   h = loglog( ScatUuuStim(bmu), ScatAuuStim(bmu), 'go'); hold on;
   set(h,'Markersize',5);
   h = loglog( ScatUuuStim(bio), ScatAuuStim(bio), 'g.'); hold on;
   set(h,'Markersize',15);  
 end
 loglog( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 
 axis tight;
 V = axis;
 h = xlabel('Away Rate (hz)');
 set(h,'FontSize',14);
 h = ylabel('Towards Rate (hz)');
 set(h,'FontSize',14);
 h = title('Stim-Locked Window');
 set(h,'FontSize',14);
 h = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.8),sprintf('N=%d',N));
 set(h,'FontSize',14);
 %**** test for significant difference 
 pval = signrank(ScatUuuStim,ScatAuuStim);
 h = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.6),sprintf('p=%6.4f',pval));
 set(h,'FontSize',14);
 %*****************************************
 
 %***** scatter plot stim-locked
 subplot('position',[0.6 0.2 0.3 0.6]);
 plot(ScatAuuPre,ScatUuuPre,'k.'); 
 maxo = max(max(ScatAuuPre),max(ScatUuuPre)) * 1.2;
 mino = min(min(ScatAuuPre),min(ScatUuuPre)) * 0.8;
 if (1)
   h = loglog( ScatUuuPre, ScatAuuPre, 'k.'); hold on;
   set(h,'Markersize',18);
 else
   h = loglog( ScatUuuPre(amu), ScatAuuPre(amu), 'ko'); hold on;
   set(h,'Markersize',5);
   h = loglog( ScatUuuPre(aio), ScatAuuPre(aio), 'k.'); hold on;
   set(h,'Markersize',18);
   h = loglog( ScatUuuPre(bmu), ScatAuuPre(bmu), 'go'); hold on;
   set(h,'Markersize',5);
   h = loglog( ScatUuuPre(bio), ScatAuuPre(bio), 'g.'); hold on;
   set(h,'Markersize',18);    
 end
 loglog( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 axis tight;
 V = axis;
 h = xlabel('Away Rate (hz)');
 set(h,'FontSize',14);
 h = ylabel('Towards Rate (hz)');
 set(h,'FontSize',14);
 h = title('Sac-Locked Window');
 set(h,'FontSize',14);
 h = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.8),sprintf('N=%d',N));
 set(h,'FontSize',14);
 %**** test for significant difference 
 pval = signrank(ScatUuuPre,ScatAuuPre);
 h = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.6),sprintf('p=%6.4f',pval));
 set(h,'FontSize',14);