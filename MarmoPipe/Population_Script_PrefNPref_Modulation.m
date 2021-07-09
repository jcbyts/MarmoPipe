%**** POPULATION SCRIPT 
%****** runs analysis across the whole population
%****** and any stats needed

%% ***** get files from Info folder *****
athome = 0;
if (athome == 1)
   INFO_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\InfoSU_JM';
else
   %INFO_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\Info';
   %INFO_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\Info4';
   INFO_DATA_DIR = 'C:\Users\Shanna\Dropbox\FovealTrans\DataAnalysis\Info4';
end

%********* Define some analysis Windows
StimWin = [50,150];
PreWin = [-50,0];
DSI_THRESH = 0.06;
%*** setup structs to store read-in data
PopAuuPrePref = [];
PopUuuPrePref = [];
PopAuuPreNPref = [];
PopUuuPreNPref = [];
%**********
NPopAuuPrePref = [];
NPopUuuPrePref = [];
NPopAuuPreNPref = [];
NPopUuuPreNPref = [];
%******************
%*** loop over Info files and collect desired data
xdir = dir(INFO_DATA_DIR);
N = 0;
for k = 3:size(xdir,1)
       filename = [INFO_DATA_DIR,filesep,xdir(k).name];
       disp(sprintf('Reading info file %s',filename));
       load(filename);
       %***** criteria for inclusion
       if (Info.DSI <= DSI_THRESH)
           continue;
       end
       %**********   better, just removed file from database (too few away
       %**********       trials to get any realistic estimates)
%        if (k == 6)
%            disp('***********');
%            disp('temp hack, removed one bad file, number 6');
%            disp(filename);
%            disp('***********');
%            continue; % skip unit that crashes
%        end
       
       PopAuuPrePref = [PopAuuPrePref ; Info.AuuPrePref];
       PopUuuPrePref = [PopUuuPrePref ; Info.UuuPrePref];
       PopAuuPreNPref =  [PopAuuPreNPref  ; Info.AuuPreNPref];
       PopUuuPreNPref =  [PopUuuPreNPref  ; Info.UuuPreNPref];
       %*******
       MaxA = max( [Info.AuuPrePref Info.AuuPreNPref]);
       MaxU = max( [Info.UuuPrePref Info.UuuPreNPref]);
       NPopAuuPrePref = [NPopAuuPrePref ; (Info.AuuPrePref * (100/MaxA))];
       NPopUuuPrePref = [NPopUuuPrePref ; (Info.UuuPrePref * (100/MaxU))];
       NPopAuuPreNPref =  [NPopAuuPreNPref  ; (Info.AuuPreNPref * (100/MaxA))];
       NPopUuuPreNPref =  [NPopUuuPreNPref  ; (Info.UuuPreNPref * (100/MaxU))];
       %*****************************
       BefStim = Info.BefStim;
       AftStim = Info.AftStim;
       BefPre = Info.BefSac;
       AftPre = Info.AftSac;
       %*********
       N = N + 1;
       clear Info;
end

%******* compute the mean stats for average normalized PSTH
AuuPrePref = mean(NPopAuuPrePref);
AsuPrePref = std(NPopAuuPrePref)/sqrt(N);
UuuPrePref = mean(NPopUuuPrePref);
UsuPrePref = std(NPopUuuPrePref)/sqrt(N);
AuuPreNPref = mean(NPopAuuPreNPref);
AsuPreNPref = std(NPopAuuPreNPref)/sqrt(N);
UuuPreNPref = mean(NPopUuuPreNPref);
UsuPreNPref = std(NPopUuuPreNPref)/sqrt(N);

 %******* plot the results now
 H = figure;
 set(H,'Position',[100 100 1200 600]);
 %****** general plot params
 SM = 2;
 Vmax = max(max(AuuPrePref+(SM*AsuPrePref)),max(UuuPrePref+(SM*UsuPrePref)));
 Vmax = Vmax * 1.1;
 
 %******* plot stim locked PSTH
 subplot('position',[0.1 0.2 0.35 0.6]);
 %****** subset of trials where no ori change occured
 h2 = plot(-BefPre:AftPre,AuuPrePref,'g-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(-BefPre:AftPre,AuuPrePref+(SM*AsuPrePref),'g-');hold on;
 h2b = plot(-BefPre:AftPre,AuuPrePref-(SM*AsuPrePref),'g-');hold on;
 %****** subset of trials where ori change DID occur
 h2 = plot(-BefPre:AftPre,AuuPreNPref,'b-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(-BefPre:AftPre,AuuPreNPref+(SM*AsuPreNPref),'b-');hold on;
 h2b = plot(-BefPre:AftPre,AuuPreNPref-(SM*AsuPreNPref),'b-');hold on;
 %************
 axis tight;
 ax = gca;
 c = ax.FontSize;
 ax.FontSize = 14;
 V = axis;
 %axis([-BefPre AftPre 0 Vmax]);
 axis([-250 AftPre 0 Vmax]);
 plot([0,0],[0,Vmax],'k-');
 plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
 plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
 t = xlabel('Time (ms)');
 t.FontSize = 18;
 t = ylabel('Norm. Rate (hz)');
 t.FontSize = 18;
 t = title('Towards');
 t.FontSize = 22;
 t.Color = 'red';
 %t = text((V(1)+(V(2)-V(1))*0.2),(V(3)+(V(4)-V(3))*0.7),sprintf('N=%d',N));
 t = text(-225,80,sprintf('N=%d',N));
 t.FontSize = 13;
 %******************
 
 %******* plot stim locked PSTH
 subplot('position',[0.55 0.2 0.35 0.6]);
 %****** subset of trials where no ori change occured
 %****** subset of trials where no ori change occured
 h2 = plot(-BefPre:AftPre,UuuPrePref,'g-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(-BefPre:AftPre,UuuPrePref+(SM*UsuPrePref),'g-');hold on;
 h2b = plot(-BefPre:AftPre,UuuPrePref-(SM*UsuPrePref),'g-');hold on;
 %****** subset of trials where ori change DID occur
 h2 = plot(-BefPre:AftPre,UuuPreNPref,'b-');hold on;
 set(h2,'Linewidth',2);
 h2b = plot(-BefPre:AftPre,UuuPreNPref+(SM*UsuPreNPref),'b-');hold on;
 h2b = plot(-BefPre:AftPre,UuuPreNPref-(SM*UsuPreNPref),'b-');hold on;
 %************
 axis tight;
 ax = gca;
 c = ax.FontSize;
 ax.FontSize = 14;
 V = axis;
 %axis([-BefPre AftPre 0 Vmax]);
 axis([-250 AftPre 0 Vmax]);
 plot([0,0],[0,Vmax],'k-');
 plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
 plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
 t = xlabel('Time (ms)');
 t.FontSize = 18;
 t = ylabel('Norm. Rate (hz)');
 t.FontSize = 18;
 t = title('Away');
 t.FontSize = 22;
 %t = text((V(1)+(V(2)-V(1))*0.2),(V(3)+(V(4)-V(3))*0.7),sprintf('N=%d',N));
 t = text(-225,80,sprintf('N=%d',N));
 t.FontSize = 13;

 