function h = initFig(varargin)
h = figure('Color','white');

if nargin>=1
    stdSize = varargin{1};
else
    stdSize = [13.2 22.23];
end
% Plos compBio: Maximum width and height [19.05 (13.2), 22.50]

Args = struct('Color','white','Units','centimeters',...
              'PapSize',stdSize ,'PapPos', [0 0 stdSize],...
              'FigPos', [0 0 stdSize],'Orientation','landscape');
if nargin>=2
          Args = myParseArgs(varargin{2},  Args);
end
% Args = parseArgs(inArgs, Args); ... % fill the arg-struct with values entered by the user
% %           {'Holdaxis'}, ... %this argument has no value (flag-type)
% %           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});

% Args.PapSize
set(h,'PaperOrientation',Args.Orientation);
set(h,'PaperUnits',Args.Units);
set(h,'PaperSize', Args.PapSize);
set(h,'PaperPosition', Args.PapPos);
set(h,'Units',Args.Units);
set(h,'OuterPosition',Args.FigPos);

end
