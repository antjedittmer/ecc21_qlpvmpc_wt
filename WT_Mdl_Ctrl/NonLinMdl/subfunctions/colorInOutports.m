function colorInOutports(aSimulinkMdl)
% color in and outports of a model

if ~nargin
    aSimulinkMdl = bdroot;
    if isempty(aSimulinkMdl)
        warning('No Simulink model open')
    end
end

% BlockType, 'Inport'; BackgroundColor, 'magenta'
aBlockType = 'Inport';
aBackgroundColor =  'magenta';
colorInports(aSimulinkMdl, aBlockType, aBackgroundColor)

% BlockType, 'Outport'; BackgroundColor, 'cyan'
aBlockType = 'Outport';
aBackgroundColor =  'cyan';
colorInports(aSimulinkMdl, aBlockType, aBackgroundColor)

function colorInports(aSimulinkMdl, aBlockType, aBackgroundColor)

gbInports = find_system(aSimulinkMdl,...           % search vdp for all Gain blocks,
    'LookUnderMasks','all','BlockType', aBlockType);
for idx = 1: length(gbInports)
    set_param(gbInports{idx},'BackgroundColor',aBackgroundColor)
end