function [strExtr,legCell] = createLegendCell(modelsL,firstCellString)
% setLegCell creates a cell with the model names for the legend. 
% It also creates a colored title legend.

if nargin < 2 || isempty(firstCellString)
    firstCellString = {'Bianchi 8 states'};
end

legCell = [firstCellString;modelsL];
cl = colormap('lines');
strExtr = '';
for ilegC = 1: length(modelsL)
    strExtr = [strExtr,...
        ', \color[rgb]{',num2str(cl(ilegC,:)),'}',modelsL{ilegC}]; %#ok<AGROW> This fast
end