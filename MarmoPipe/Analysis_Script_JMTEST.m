%*********** ANALYSIS SCRIPT *******************************
%***********************************************************
%*** NOTE:  The goal here is just to analyze a single file
%***        Reads in a Exp struct created from the ImportScript
%***        that would be placed in the data folder, Exp.mat
%************************************************************

%% Spike File importing
%******* Select the folder for analysis *********************
%FileTag = 'Ellie_301118';
% FileTag = 'Ellie_291118';
% FileTag = 'Ellie_041218';
% FileTag = 'Ellie_211118';
% FileTag = 'Ellie_081218';
% FileTag = 'Ellie_121218';
%FileTag = 'Ellie_151218';
%FileTag = 'Ellie_161218';
%FileTag = 'Ellie_171218';
FileTag = 'Ellie_251218';  % to build spike sort
%FileTag = 'Ellie_261218';
%FileTag = 'Ellie_211218';
%FileTag = 'Ellie_080119';
%FileTag = 'Ellie_090119';
%FileTag = 'Ellie_140119';
%FileTag = 'Ellie_150119';
%FileTag = 'Ellie_170119';
%FileTag = 'Ellie_190119';
%FileTag = 'Ellie_200119';
%FileTag = 'Ellie_210119';
%FileTag = 'Ellie_230119';
%FileTag = 'Ellie_240119';

%FileTag = 'Ellie_290119';
%FileTag = 'Ellie_300119';
%FileTag = 'Ellie_310119';
%FileTag = 'Ellie_010219';
%FileTag = 'Ellie_050219';
%FileTag = 'Ellie_060219';
%FileTag = 'Ellie_070219';
%FileTag = 'Ellie_080219';
%FileTag = 'Ellie_090219';
%FileTag = 'Ellie_100219';
%FileTag = 'Ellie_110219';
%FileTag = 'Ellie_120219';
%FileTag = 'Ellie_150219d';
%FileTag = 'Ellie_150219e';
%FileTag = 'Ellie_160219';
%FileTag = 'Ellie_160219b';
%FileTag = 'Ellie_170219';
%FileTag = 'Ellie_170219b';
%FileTag = 'Ellie_180219';
%FileTag = 'Ellie_180219b';
%FileTag = 'Ellie_260219';
%FileTag = 'Ellie_270219b';
%FileTag = 'Ellie_270219c';
%FileTag = 'Ellie_020319b';
%FileTag = 'Ellie_030319';
%FileTag = 'Ellie_040319';
%FileTag = 'Ellie_040319c';
%FileTag = 'Ellie_050319d';
%FileTag = 'Ellie_050319c';
%FileTag = 'Ellie_060319';
%FileTag = 'Ellie_070319b';
%FileTag = 'Ellie_080319b';
%FileTag = 'Ellie_080319c';
%FileTag = 'Ellie_110319';
%FileTag = 'Ellie_110319b';
%FileTag = 'Ellie_120319c';
%FileTag = 'Ellie_120319d';
%FileTag = 'Ellie_140319b';
%FileTag = 'Ellie_140319c';
%FileTag = 'Ellie_210319';
%FileTag = 'Ellie_220319';
%FileTag = 'Ellie_230319'; % Foveal RF
%FileTag = 'Ellie_240319'; 
% FileTag = 'Ellie_250319'; 
%FileTag = 'Ellie_260319'; 
%FileTag = 'Ellie_020419';  % foveal RF, mo map too
%FileTag = 'Ellie_030419'; 
%FileTag = 'Ellie_040419'; % No RF
%FileTag = 'Ellie_050419';
%FileTag = 'Ellie_060419';
%FileTag = 'Ellie_080419';
%FileTag = 'Ellie_130419';
%FileTag = 'Ellie_230419';
%FileTag = 'Ellie_240419b';
%FileTag = 'Ellie_260419';
%FileTag = 'Ellie_260419b';
% FileTag = 'Ellie_290419b';
SPClust = 1;
% SPClust = 3;
%************* 
athome = 0;
if (athome == 1)
   PROCESSED_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\Processed';
else
   %PROCESSED_DATA_DIR = 'C:\Users\mitchelllab\Dropbox\FovealTrans\DataAnalysis\Processed';
   PROCESSED_DATA_DIR = 'C:\Users\Shanna\Dropbox\FovealTrans\DataAnalysis\Processed';
end

%% Grab Exp File that was previously imported
ExpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'.mat'];
disp(sprintf('Loading processed Exp file, %s',FileTag));
load(ExpFile);  % struct Exp will be loaded
disp('File loaded');

%% Analyze the spatial RF
SPClust = 1; 
if (0)
    if (1)
      GRID.box = [0,0,20,20];  % center x,y, width and height
      GRID.div = 2.0;          % size of divisions
    else
      GRID.box = [0,0,40,40];  % center x,y, width and height
      GRID.div = 4.0;          % size of divisions       
    end
else
  GRID.box = [0,0,14,14];  % center x,y, width and height
  GRID.div = 1.0;
end
[StimX,StimXP,StimXN,StimY] = Forage.StimMatrix_ForageSpatialKernel(Exp,GRID,SPClust);  %step 1 
Forage.PlotForageSpatialKernel(StimX,StimY,GRID,FileTag);     % step 2, R
% Forage.PlotForageSpatialKernel(StimXP,StimY,GRID,FileTag);  
% Forage.PlotForageSpatialKernel(StimXN,StimY,GRID,FileTag);  

%% Motion Analysis of Space RF
SPClust = 1; 
GRID.box = [0,0,40,40];  % center x,y, width an d height
GRID.div = 4.0; 
[MoStimX,MoStimY] = Forage.StimMatrix_ForageMotionSpatialKernel(Exp,GRID,SPClust);  %step 1 
Forage.PlotForageMotionSpatialKernel(MoStimX,MoStimY,GRID,FileTag);     % step 2, R

%% Analyze the foraging experiment to label saccade onsets in same frame timing
% RUN THIS FOR FORAGE4
[SacX] = Forage.SacMatrix_ForageGratingKernel(Exp);

%% User Defined Analysis to examine reverse correlation kernel
%*****************************************************************
% looks at temporal kernel to gratings, plus spat freq x orientation
%*******
% RUN THIS FOR FORAGE4 for each unit, seperated steps, but all most be
% rerun per unit (steps 2 or 3 can be run after step 1)
SPClust = 1;
[StimX,StimY] = Forage.StimMatrix_ForageGratingKernel(Exp,SPClust);  %step 1
Forage.PlotForageGratingKernel(StimX,StimY,FileTag);     % step 2, RF

%% if you want to looks at saccade triggered rate average and tuning
% last param, 1 - sac onset aligned, 4- sac offset, 5-stim reappear fovea
Forage.PlotForageSaccadeKernel(SacX,StimX,StimY,FileTag,1);  % step 3, sac trig

%%  Run plot of motion direction selectivity
%SPClust = 1;
Forage.PlotForageMotionSelectivity(Exp,SPClust);  % plot raster by motion direction

%% User Defined Analysis to examine reverse correlation kernel of CamoFlag
%*****************************************************************
% first pass, would just like to look at temporal response to
% a flash from an oriented grating, and maybe orientation selectivity
[StimX,StimY] = Flag.StimMatrix_FlagGratingKernel(Exp,SPClust);
Flag.PlotFlagGratingKernel(StimX,StimY,FileTag);

%% User Defined Analysis to examine reverse correlation kernel from  
% the CamoFlag experiments that have background Hartley stimuli
%*****************************************************************
% looks at temporal kernel to gratings, plus spat freq x orientation
[StimX,StimY] = Flag.StimMatrix_CamoGratingKernel(Exp);
Forage.PlotForageGratingKernel(StimX,StimY,FileTag);

%% User Defined Analysis to examine saccades to targets and neural activity
%*****************************************************************
% Flag.PlotActivityAroundSaccade(Exp,TimeLock,ShowRaw)
%    Exp - the data struct
%    Time-lock : 0 - saccade onset, 1 - saccade offset  2 - stim onset
%    ShowRaw - if 1, show the raw counts instead of mean and error bars
SPClust = 1;
Flag.PlotActivityAroundSaccade(Exp,1,0,SPClust);

%% User Defined Analysis to examine saccades to targets and neural activity
%*****************************************************************
% Flag.PlotActivityAroundSaccade(Exp,TimeLock,ShowRaw)
%    Exp - the data struct
%    Time-lock : 0 - saccade onset, 1 - saccade offset  2 - stim onset
%    ShowRaw - if 1, show the raw counts instead of mean and error bars
SPClust = 1;
Flag.PlotMoActivityAroundSaccade(Exp,2,0,SPClust);

%% User Defined Analysis to examine saccades to targets and neural activity
%*****************************************************************
% Flag.PlotActivityAroundSaccade(Exp,TimeLock,ShowRaw)
%    Exp - the data struct
SPClust = 1;
Flag.PlotMoBeforeSaccade(Exp,SPClust);
%***********************

%% User Defined Analysis to examine saccades to targets and neural activity
%*****************************************************************
SPClust = 1;
GRID.box = [0,0,20,20];  % center x,y, width and height
GRID.div = 2.0;          % size of divisions
if (0)
   [StimX,StimXP,StimXN,StimY] = Forage.StimMatrix_ForageSpatialKernel(Exp,GRID,SPClust);
else
   [StimX,StimY] = Forage.StimMatrix_ForageMotionSpatialKernel(Exp,GRID,SPClust);  %step 1 
end

% Forage.PlotForageSpatialKernel(StimX,StimY,GRID,FileTag);     % step 2, R
% Forage.PlotForageMotionSpatialKernel(StimX,StimY,GRID,FileTag);     % step 2, R

%%
Info = Flag.PlotMoBeforeSaccade_RF(Exp,SPClust,StimX,StimY,GRID,FileTag);
%%
Flag.PlotInfo_MoBeforeSaccade_PrefNPref_ROC(Info);
%%
Flag.PlotInfo_MoBeforeSaccade_VonMises_Fit(Info);
%%
Flag.PlotInfo_MoBeforeSaccade_Rate_Raster(Info);
%%
%*******************
InfoFile = [PROCESSED_DATA_DIR,filesep,FileTag,sprintf('%d',SPClust),'_I.mat'];
save(InfoFile,'Info');
%******* save the file
h = gcf;
z = getframe(h);
FigFile = [PROCESSED_DATA_DIR,filesep,FileTag,sprintf('%d',SPClust),'_F.tiff'];
imwrite(z.cdata,FigFile,'Compression','none');
% %***********************


