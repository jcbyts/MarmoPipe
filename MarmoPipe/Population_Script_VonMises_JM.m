%**** POPULATION SCRIPT 
%****** runs analysis across the whole population
%****** and any stats needed

% ***** get files from Info folder *****
athome = 0;
if (athome == 1)
   INFO_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\Info2'; %Wid_JM';
else
   %INFO_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\InfoWid_JM';
   %INFO_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\Info2';
   INFO_DATA_DIR = 'C:\Users\Shanna\Dropbox\FovealTrans\DataAnalysis\Info4';
end

%********* Define some analysis Windows
DSI_THRESH = 0.1;
FVAL_THRESH = -10e2; % Goodness of fit
ISO_THRESH = 0;  % zero if include all
WAVE_THRESH = 0; % zero if include ll
%*** setup structs to store read-in data
Abase = [];
Ubase = [];
Aamp = [];
Uamp = [];
Awid = [];
Uwid = [];
DSI = [];
Adsi = [];
Udsi = [];
%*** loop over Info files and collect desired data
xdir = dir(INFO_DATA_DIR);
N = 0;
tiny = 0.00001;
 for k = 3:size(xdir,1)
       filename = [INFO_DATA_DIR,filesep,xdir(k).name];
       load(filename);
       %****** waveform
       wave = Info.Wave;
       minwave = min(wave);
       za = find( wave == minwave);
       maxwave = max(wave(za(1):end));
       zb = find( wave == maxwave);
       wdiff = zb(1) - za(1);
       %**********
       if (WAVE_THRESH ~= 0)
         if (WAVE_THRESH > 0)
           if (wdiff < WAVE_THRESH)
            continue;
           end
         else
           if (wdiff >= -WAVE_THRESH)
            continue;
           end    
         end
       end
       %***** criteria for inclusion
   
       if (Info.DSI <= DSI_THRESH)
           continue;
       end
       if (Info.Tfvalstim > FVAL_THRESH)
           continue;
       end
       if (Info.Iso < ISO_THRESH)
           continue;
       end
       
       %*******
       % figure(10); hold off;
       % plot(wave,'b-'); hold on;
       % wdiff
       %******
       % wdiff
       % input('check');
       %*********
       disp(sprintf('Reading sample %d info file %s',(N+1),xdir(k).name));
       DSI = [DSI ; Info.DSI];
       %*******
       if (0)
           Abase = [Abase ; [Info.Afitstim.mu(1),Info.Afitstim.sem(1)]];
           Ubase = [Ubase ; [Info.Ufitstim.mu(1),Info.Ufitstim.sem(1)]];
           Aamp = [Aamp ; [Info.Afitstim.mu(2),Info.Afitstim.sem(2)]];
           Uamp = [Uamp ; [Info.Ufitstim.mu(2),Info.Ufitstim.sem(2)]];
           Awid = [Awid ; [Info.Afitstim.mu(3),Info.Afitstim.sem(3)]];
           Uwid = [Uwid ; [Info.Ufitstim.mu(3),Info.Ufitstim.sem(3)]];
           Adsi = [Adsi ; [Info.ADSI]];
           Udsi = [Udsi ; [Info.UDSI]];
       else
           Abase = [Abase ; [Info.Afitpre.mu(1),Info.Afitpre.sem(1)]];
           Ubase = [Ubase ; [Info.Ufitpre.mu(1),Info.Ufitpre.sem(1)]];
           Aamp = [Aamp ; [Info.Afitpre.mu(2),Info.Afitpre.sem(2)]];
           Uamp = [Uamp ; [Info.Ufitpre.mu(2),Info.Ufitpre.sem(2)]];
           Awid = [Awid ; [Info.Awidpre.mu,Info.Awidpre.sem]];
           Uwid = [Uwid ; [Info.Uwidpre.mu,Info.Uwidpre.sem]];
           Adsi = [Adsi ; [Info.ADSI]];
           Udsi = [Udsi ; [Info.UDSI]];
%           Awid = [Awid ; [Info.Afitpre.mu(3),Info.Afitpre.sem(3)]];
%           Uwid = [Uwid ; [Info.Ufitpre.mu(3),Info.Ufitpre.sem(3)]];
       end
       %*********
       N = N + 1;
       clear Info;
 end

 %******* plot the results now
 H = figure;
 set(H,'Position',[100 100 800 800]);

 %***** now perform a scatter plot *********
 %***** scatter plot stim-locked
 subplot('position',[0.1 0.1 0.35 0.35]);
 %h = plot(Ubase(:,1),Abase(:,1),'k.'); hold on;
 %set(h,'Markersize',15);
 MaxBaseVar = max(max(Abase(:,2)),max(Ubase(:,2)));
 for k = 1:size(Ubase,1)
     bv = 0.9 * (0.5*(Abase(k,2)+Ubase(k,2)))/MaxBaseVar;
     colo = [bv,bv,bv];
     h = plot(Ubase(k,1),Abase(k,1),'k.'); hold on;
     set(h,'Color',colo);
     set(h,'Markersize',15);
     h = plot([Ubase(k,1),Ubase(k,1)],[(Abase(k,1)-Abase(k,2)),(Abase(k,1)+Abase(k,2))],'k-');
     set(h,'Color',colo);
     h = plot([(Ubase(k,1)-Ubase(k,2)),(Ubase(k,1)+Ubase(k,2))],[Abase(k,1),Abase(k,1)],'k-');
     set(h,'Color',colo);
 end
 maxo = max(max(Abase),max(Ubase)) * 1.1;
 mino = min(min(Abase),min(Ubase)) * 0.2;
 plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
%labels = cellstr(num2str([1:N]'));
 axis tight;
 V = axis;
 h = xlabel('Away');
 set(h,'FontSize',20);
 h = ylabel('Towards');
 set(h,'FontSize',20);
 set(h,'Color',[1,0,0]);
 t = title('Baseline');
 t.FontSize = 16;
 t = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.90),sprintf('N=%d',N));
 t.FontSize = 12;
 %**** test for significant difference 
 pval = signrank(Ubase(:,1),Abase(:,1));
 t = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.80),sprintf('p=%6.4f',pval));
 t.FontSize = 12;
 
 %************************
 subplot('position',[0.6 0.1 0.35 0.35]);
 %h = plot(Uamp(:,1),Aamp(:,1),'k.'); hold on;
 %set(h,'Markersize',15);
 MaxAmpVar = max(max(Aamp(:,2)),max(Uamp(:,2)));
 for k = 1:size(Ubase,1)
     bv = 0.9 * (0.5*(Aamp(k,2)+Uamp(k,2)))/MaxAmpVar;
     colo = [bv,bv,bv];
     h = plot(Uamp(k,1),Aamp(k,1),'k.'); hold on;
     set(h,'Color',colo);
     set(h,'Markersize',15);
     h = plot([Uamp(k,1),Uamp(k,1)],[(Aamp(k,1)-Aamp(k,2)),(Aamp(k,1)+Aamp(k,2))],'k-');
     set(h,'Color',colo);
     h = plot([(Uamp(k,1)-Uamp(k,2)),(Uamp(k,1)+Uamp(k,2))],[Aamp(k,1),Aamp(k,1)],'k-');
     set(h,'Color',colo);
 end
 maxo = max(max(Aamp),max(Uamp)) * 1.1;
 mino = min(min(Aamp),min(Uamp)) * 0.2;
 plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 %************
 axis tight;
 V = axis;
 h = xlabel('Away');
 set(h,'FontSize',20);
 h = ylabel('Towards');
 set(h,'Color',[1,0,0]);
 set(h,'FontSize',20);
 t = title('Amplitude');
 t.FontSize = 16;
 t = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.90),sprintf('N=%d',N));
 t.FontSize = 12;
 %**** test for significant difference 
 pval = signrank(Uamp(:,1),Aamp(:,1));
 t = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.80),sprintf('p=%6.4f',pval));
 t.FontSize = 12;

 
 %************************
 subplot('position',[0.1 0.6 0.35 0.35]);
 %h = plot(log(Uwid),log(Awid),'k.'); hold on;
 %set(h,'Markersize',15);
 MaxWidVar = max(max(Awid(:,2)),max(Uwid(:,2)));
 for k = 1:size(Uwid,1)
    bv = 0.9 * (0.5*(Awid(k,2)+Uwid(k,2)))/MaxWidVar;
     colo = [bv,bv,bv];
     h = plot(Uwid(k,1),Awid(k,1),'k.'); hold on;
     set(h,'Color',colo);
     set(h,'Markersize',15);
     h = plot([Uwid(k,1),Uwid(k,1)],[(Awid(k,1)-Awid(k,2)),(Awid(k,1)+Awid(k,2))],'k-');
     set(h,'Color',colo);
     h = plot([(Uwid(k,1)-Uwid(k,2)),(Uwid(k,1)+Uwid(k,2))],[Awid(k,1),Awid(k,1)],'k-');
     set(h,'Color',colo);
 end
 axis tight;
 V = axis;
 maxo = 360; %max(V);
 mino = 0; %min(V);
 plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 axis([mino maxo mino maxo]);
 %**********
 %axis tight;
 V = axis;
 h = xlabel('Away');
 set(h,'FontSize',20);
 h = ylabel('Towards');
 set(h,'Color',[1,0,0]);
 set(h,'FontSize',20);
 t = title('Half-Width');
 t.FontSize = 16;
 t = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.90),sprintf('N=%d',N));
 t.FontSize = 12;
 %**** test for significant difference 
 pval = signrank(Uwid(:,1),Awid(:,1));
 t = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.80),sprintf('p=%6.4f',pval));
 t.FontSize = 12;
 
 
 %************************
 subplot('position',[0.6 0.6 0.35 0.35]);
 %h = plot(log(Uwid),log(Awid),'k.'); hold on;
 %set(h,'Markersize',15);
 for k = 1:size(Udsi,1)
     h = plot(Udsi(k,1),Adsi(k,1),'k.'); hold on;
     set(h,'Color',colo);
     set(h,'Markersize',15);
 end
 axis tight;
 V = axis;
 maxo = max(V) * 1.25;
 mino = min(V) * 0.75;
 plot( [mino,maxo], [mino,maxo], 'k--');  % line of unity
 %**********
 axis tight;
 V = axis;
 h = xlabel('Away');
 set(h,'FontSize',20);
 h = ylabel('Towards');
 set(h,'Color',[1,0,0]);
 set(h,'FontSize',20);
 t = title('DSI');
 t.FontSize = 16;
 t = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.90),sprintf('N=%d',N));
 t.FontSize = 12;
 %**** test for significant difference 
 pval = signrank(Udsi(:,1),Adsi(:,1));
 t = text((V(1)+(V(2)-V(1))*0.1),(V(3)+(V(4)-V(3))*0.80),sprintf('p=%6.4f',pval));
 t.FontSize = 12;
 
 
 
 