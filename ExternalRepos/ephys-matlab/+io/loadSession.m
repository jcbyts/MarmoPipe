function [sessionInfo, ops, info] = loadSession(oepath)
% LOADSESSION loads the meta data about an OE/PLDAPS session
% Inputs:
%   oepath@string    - path to session
%   oepath can also be a meta data entry
% Outputs:
% 	sessionInfo@struct - info about session(s)
%   ops@struct         - ops struct. built by io.oe2dat --> used for KiloSort
%   info@struct        - info struct. built by io.oe2dat --> meta-data required to make sense of raw ephys samples
% Example Call:
%   oepath = 'C:\Data\Ellie_2017-08-09_13-04-23_ShankD15MT6';
%   [sessionInfo, ops, info] = loadSession(oepath)
%
%   meta = io.getExperimentsAnd('Subject', 'Ellie', 'StimulusProtocols', 'hartleyFF');
%   [sessionInfo, ops, info] = io.loadSession(meta(1,:))

SERVER_DATA_DIR = 'Z:\Data\PLDAPS\Ellie';

if nargin == 0
    oepath = uigetdir(SERVER_DATA_DIR); %getpref('EPHYS', 'SERVER_DATA'));
end

if istable(oepath)
    oepath = fullfile(SERVER_DATA_DIR, oepath.Directory{1});
     % getpref('EPHYS', 'SERVER_DATA'), oepath.Directory{1});
end

ops = io.loadOps(oepath);
info = io.loadEphysInfo(oepath);

load(fullfile(oepath, 'session_info.mat'))