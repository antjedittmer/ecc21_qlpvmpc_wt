function colorInOutports(aSimulinkMdl,aBackgroundColorOut,aBackgroundColorIn)
% colorInOutports colors in and outports of a model
%
% Inputs: * aSimulinkMdl: model to be changed. If omitted, the current
%           model is used
%         * aBackgroundColorOut: Outport color. If omitted, set to 'magenta'
%         * aBackgroundColorIn: Inport color. If omitted, set to 'cyan'

% Antje Dittmer, DLR 

%% Set default model, out-color, in-color
if ~nargin
    aSimulinkMdl = bdroot;
    if isempty(aSimulinkMdl)
        warning('No Simulink model open. Return without action')
        return;
    end
end

if nargin < 2 || isempty(aBackgroundColorOut)
    aBackgroundColorOut =  'magenta';
end

if nargin < 3 || isempty(aBackgroundColorIn)
    aBackgroundColorIn =  'cyan';
end

%% BlockType, 'Inport'; BackgroundColor, 'magenta'
aBlockType = 'Inport';
colorInports(aSimulinkMdl, aBlockType, aBackgroundColorOut)

%% BlockType, 'Outport'; BackgroundColor, 'cyan'
aBlockType = 'Outport';

colorInports(aSimulinkMdl, aBlockType, aBackgroundColorIn)

function colorInports(aSimulinkMdl, aBlockType, aBackgroundColor)

gbInports = find_system(aSimulinkMdl,...           % search vdp for all Gain blocks,
    'LookUnderMasks','all','BlockType', aBlockType);
for idx = 1: length(gbInports)
    set_param(gbInports{idx},'BackgroundColor',aBackgroundColor)
end