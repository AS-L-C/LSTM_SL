clc
close all
clear all

%% Initialize figure-------------------------------------------------------
figH = initFig([13.2, 22.50]); %It works!
% figH = initFig([19.05, 22.50]); %It works!
%% Create layout-----------------------------------------------------------
% Set general layout
pn = struct('r',7,'c',4,'rtop',4,'rbot',3);
handles.main    = tiledlayout(figH, pn.r, 1, 'padding', 'normal','TileSpacing','none'); % Create a layout of 5 rows and 1 column
handles.main.Units = 'centimeters';
handles.main.OuterPosition = [.1 0 13.2 22.50*.9];
% Create 'handles.top' sublayout
handles.top     = tiledlayout(handles.main, pn.rtop, pn.c, 'padding', 'none','TileSpacing','none'); % Create a 4X3 sublayout 'handles.top'
handles.top.Layout.Tile     = [1]; % Make 'handles.top' start from tile 1 of 'haldles.main'
handles.top.Layout.TileSpan = [pn.rtop 1]; % Make 'handles.top' span 4 rows and 1 column
% Create 'handles.bottom' sublayout
handles.bottom  = tiledlayout(handles.main, pn.rbot, pn.c, 'padding', 'none','TileSpacing','none');  % 3X1
% handles.bottom.Layout.Tile     = [pn.rtop*pn.c + 1];
handles.bottom.Layout.Tile     = [pn.rtop + 1]; %Starts from 5th tile of handles.main
handles.bottom.Layout.TileSpan = [pn.rbot 1];
% nexttile(handles.bottom)
% plot(rand(1,10))

%% Initialize model parameters
% p = mfilename('fullpath'); %Full path of current file
load('./Manuscript Figures/Data/Parameters/s82 learningTrace/allv.mat',...
    'systpEst','estimPar','usts');
indValStates = [1 2 5 6 9 10]; %Valid states
sFuns        = struct('syspec',@syspec_ssim_v15,'ffun',@ffun_ssim_v15,'gfun',@gfun_ssim_v13);

% cols.default = get(groot,'defaultAxesColorOrder');
n.states = length(indValStates);

%% Load and simulate-------------------------------------------------------
% Load data
SData        = load('.\Manuscript Figures\Data\Experimental\IS Experiment\IS_preprocessed_v3.mat');
data         = structfun(@(x) x([2 1]),  SData.summary, 'UniformOutput', false);
clear Sdata;
n.totStrides = length(data.m{1,1});  % Same as parad{1,1}.n.ttot;

names.groups = {'SAV','INT'};
n.groups     = length(names.groups);

% Load and simulate SAV and INT paradigms----------------------------------
d              = abs(usts.extra.dest.IS);
c2 = cell(1,n.groups); parad = c2; epNames = c2; simp = c2; cout = c2;
for gr = 1:n.groups
    [parad{gr}, epNames{gr}] = loadParadigm(names.groups{gr},d);
    simp{gr} = struct('oneStep', false, 'nt', parad{gr}.n.ttot,'u', parad{gr}.u,...
        't', parad{gr}.t,  'udt2', [], 'Ts', []);
    
    % Simulate system
    cout{gr}      = mySystSim(sFuns, systpEst, simp{gr});
    parad{gr}.y   = cout{gr}.y;
    parad{gr}.z   = parad{gr}.u - parad{gr}.y;
end


%% Set other plotting parameters-------------------------------------------
fs.m   = 18;
fs.s   = 16;
fs.xs  = 14;
fs.xxs = 12;

lw.mark  = 2;
lw.lines = 3;
my_set_default(fs.xxs,lw.lines,2);
[cols.default, cols.cbr, cols.grays] = getColors();
cols.groups = [cols.default.blue; cols.default.orange];
% cols.conditions = [[1 1 1];cols.default.green + [0 .2 0];...
%                    cols.default.orange; [1 1 1];...
%                    cols.default.green + [0 .2 0];  [1 1 1]];
cols.st     = [.4 .4 .4];
cols.ea     = [.2 .2 .2];
cols.end    = [0 0 0];
cols.epochs = [cols.st; cols.ea; cols.end; cols.st; cols.end; ...
    cols.st; cols.ea; cols.end; cols.st; cols.end];
epochs.def = {[151, 155],   [201, 230],  [721,750],   [1951, 1955], ...
    [2071, 2100], [2101 2105], [2151 2180], [2671 2700], ...
    [2701 2705],  [2821, 2850]};
epochs.names = {'Start', 'Early' , 'End', ...
    'Start', 'End', ...
    'Start', 'Early' , 'End', ...
    'Start', 'End'};

%% Plot--------------------------------------------------------------------
%% A: Plot paradigms ------------------------------------------------------
sp.A = nexttile(handles.top,1,[1 pn.c]);
dels = [.035 -.035];

for gr=1:n.groups
    hold on
    prh(gr) = plot(nanTrans(parad{gr}.u + dels(gr),parad{gr}.transitions),'color', cols.groups(gr,:));
end
% Add patches
%Savings
cols.greenPatch = cols.default.green + [0 .2 0];
patchParams.sav = struct('ylims', [-d d], 'colors', repmat(cols.greenPatch, 2, 1) ,...
    'dx', 0,     'dy', .035 ,...
    'alpha', .2, 'patchInd', [2 5],...
    'patchPos', [0 0]);
addParadPatches(sp.A, parad{1}, patchParams.sav);

%Interference
patchParams.int          = patchParams.sav;
patchParams.int.colors   = cols.default.orange;
patchParams.int.patchInd = 3;
patchParams.int.patchPos = -1;
addParadPatches(sp.A, parad{2}, patchParams.int);

% Add labels
addParadLabels(sp.A,  parad{1}, fs.xxs, [-d d], [ones(1,6)], [1:6]);
addParadLabels(sp.A,  parad{2}, fs.xxs, [-d d],          -1, 3);

% Change labels
ylim([-d d])
xticks(parad{1}.transitions);
xticklabels('');
ylab(1) = ylabel('Perturb. Magn.','FontSize', fs.xxs);
yticks([-d 0 d])
% yticklabels({'2:1', '1:1', '1:2'});
yticklabels({'1', '0', '-1'});

axis tight
% legend(prh,'Savings','Interference', 'Location', 'southwest',...
%            'NumColumns',1,'Orientation', 'vertical')
% text(sp.A, 1, -d, {['\color[rgb]{' cols.default.blue '} Savings'],...
%            ['\color[rgb]{' cols.default.blue '}Interference']},'FontSize', fs.xs)
text(sp.A, 50, -d + .1 + .1 , 'Savings', 'Color', cols.default.blue, 'FontSize',   fs.xxs,'VerticalAlignment', 'bottom');
text(sp.A, 50, -d - .1 + .05 , 'Interference', 'Color', cols.default.orange, 'FontSize', fs.xxs,'VerticalAlignment', 'bottom');
% ca=gca;

% exportgraphics(figH, 'outputTest1.pdf', 'ContentType', 'vector');

%% B: time courses
plPars = struct('type','data','cols',cols,'n',n,'ylimits',[-.85 .85],...
    'patchParams',patchParams,'epochs',epochs,'d',d,'parad',{parad},...
    'fs', fs, 'expandStart', true, 'expandStartMagn', 3,...
    'addEpochLabels', true, 'rectangles', true,...
    'dotSize', 30, 'rect', '');
sp.B   = nexttile(handles.top,pn.c + 1,[2 pn.c]);
out    = plotTimeCourse(data, plPars, sp.B);



%% C: zoom in on relevant epochs
zpars.zoomIns  = {151:250,2101:2200,1951:2100-3,2701:2850-3};
zpars.zoomInsL = {'151','250';'2101','2200';'1951','2100';'2701','2850'};
zpars.zoomep   = {'A_1', 'A_2', 'N_s', 'W_o'};
zpars.legloc   = {'mr','mr','bm','mr'};
sp.C           = nexttile(handles.top, pn.c*3 + 1,[1 1]);
plPars.pn      = pn;
plPars.rect    = out.rect;
plPars.rect.fNormHeight = out.rect.normHeight;
plPars.pn.prevRows= 3;
plotZoomIN(data, plPars, zpars, handles.top)

%% D: Model time courses
plPars.type = 'model';
plPars.addEpochLabels = false;
sp.D        = nexttile(handles.bottom,1,[2 pn.c]);

plPars.rect.fNormHeight = out.rect.normHeight*(sp.D.Position(4)/sp.B.Position(4));
out                     = plotTimeCourse(data, plPars, sp.D);

%% E:
sp.E        = nexttile(handles.bottom, pn.c*2 + 1,[1 1]);
plPars.pn.prevRows= 2;
plPars.rectangles = false;
% plPars.rect.fNormHeight = out.rect.normHeight*(sp.E.Position(4)/sp.B.Position(4));
plotZoomIN(data, plPars, zpars, handles.bottom)
sp.E.XLabel.String   = 'Strides';
sp.E.XLabel.FontSize = fs.xxs';

%% Add panel labels
% Paradigm
% cy = sum(handles.top.Position([2 4])); %.9350
% cx = handles.top.Position(1) - .08;
% h  = annotation('textarrow',[cx,cx],[cy,cy],'HeadStyle', 'none', 'LineStyle', 'none','String','Paradigm','FontSize',fs.xs,'TextRotation',90);
% h  = annotation('textarrow',[cx,cx]-.01,[cy,cy]+.02,'HeadStyle', 'none', 'LineStyle', 'none','String','A','FontSize',fs.xs);
% 
% % Data
% cy = cy -.325 + .2;
% h  = annotation('textarrow',[cx,cx],[cy,cy],'HeadStyle', 'none', 'LineStyle', 'none','String','Data','FontSize',fs.xs,'TextRotation',90);
% h  = annotation('textarrow',[cx,cx]-.01,[cy,cy],'HeadStyle', 'none', 'LineStyle', 'none','String','B','FontSize',fs.xs);
% 
% % LSTM-SL
% cy = sum(handles.bottom.Position([2 4])) - .15;
% h  = annotation('textarrow',[cx,cx],[cy,cy],'HeadStyle', 'none', 'LineStyle', 'none','String','LSTM-SL','FontSize',fs.xs,'TextRotation',90);
% h  = annotation('textarrow',[cx,cx]-.01,[cy,cy]+.15,'HeadStyle', 'none', 'LineStyle', 'none','String','C','FontSize',fs.xs);
addAnnotations(figH)
function plotZoomIN(data, p, zpars, parentTile)
cols   = p.cols;
d      = p.d;
n      = p.n;
fs     = p.fs;
epochs = p.epochs;
pn     = p.pn;

n.zis  = length(zpars.zoomIns);
for zi = 1:n.zis
    cHand(zi) = nexttile(parentTile, pn.c*pn.prevRows + zi,[1 1]);
    hold on
    cs = zpars.zoomIns{zi};
    if zi==1
        delGR  = [-2 2];
    else
        delGR = [0 0];
    end

    for ig = 1:n.groups
        %             plot(cs,data.m{g}(cs));
        if strcmp(p.type,'data')
            hold on, [lh(ig), shadh(ig)] = boundedline(cs, data.m{ig}(cs), data.se{ig}(cs), ...
                '-o', 'nan', 'gap');
            shadh(ig).FaceColor    = cols.groups(ig,:);
            shadh(ig).FaceAlpha = .5;
        elseif strcmp(p.type,'model')
            if zi==1 && ig==2
                ls = '--';
            else
                ls = '-';
            end
            csla    = -p.parad{ig}.y(cs);
            csla(1) = nan;
            lh(ig)  = plot(cs, csla,'Color',...
                cols.groups(ig,:),'LineStyle',ls);
            scatter(cs(2) + delGR(ig),csla(2),p.dotSize,cols.groups(ig,:),'filled');
            
        end
        lh(ig).LineWidth = 1;
        lh(ig).Color           = cols.groups(ig,:);
        lh(ig).MarkerFaceColor = cols.groups(ig,:);
        axis tight
        
    end
    xticks([cs(1) cs(end)])
    xticklabels(zpars.zoomInsL(zi,:))
    
    if zi==3 && strcmp(p.type,'model')
        cHand(zi).YAxis.Exponent = 0;
        
        %             xtickformat('%.0f')
        cyl = cHand(zi).YLim;
        cHand(zi).YLim       = [cyl(1) .001];
        cHand(zi).YTick      = [-.003    0];
        cHand(zi).YTickLabel = {'-.003' '0'};
        %             cHand(zi).YTickLabel = ['3e-14'];
        
    end
end
% cHand(3).YTick = [.01 -.01];
% cHand(2).YLim = cHand(1).YLim;
linkaxes(cHand(1:2),'y');
% Add patches
% cHand(zi) = nexttile(handles.top, pn.c*3 + zi,[1 1]);


for zi = 1:2
    %         cHand(zi) = nexttile(handles.top, pn.c*3 + zi,[1 1]);
    cx = cHand(zi).XLim;
    cy = cHand(zi).YLim;
    hold on,   cph = patch(cHand(zi),[cx fliplr(cx)],[cy(1)*[1 1] cy(2)*[1 1]],...
        cols.greenPatch,'EdgeColor','none',...
        'FaceAlpha', p.patchParams.sav.alpha);
    uistack(cph,'bottom')
    
    
end

% Add epoch names
for zi=1:4
    cx = cHand(zi).XLim;
    cy = cHand(zi).YLim;
    myAddText(cHand(zi),cx,cy,zpars.legloc{zi},zpars.zoomep{zi},fs.xxs)
end

% Add rectangles
if p.rectangles
    rects  = {[1 2]; [1 2]; [1 3]; [1 3]};
    epZoom = {{epochs.def{1}, epochs.def{2}}; {epochs.def{6}, epochs.def{7}};...
        {epochs.def{4}, epochs.def{5}}; {epochs.def{9}, epochs.def{10}}};
    for zi = 1:4
        cx = cHand(zi).XLim;
        cy = cHand(zi).YLim;        
        %     rect.normHeight  = cax.Position(4)*rect.height/rect.range; %Height normalized with respect to figure
        %     cheight = p.rect.normHeight*diff(cy)/cHand(zi).Position(4);
        cheight = p.rect.fNormHeight*diff(cy)/cHand(zi).Position(4);
        
        ctypes = rects{zi};
        ceps   = epZoom{zi};
        for ir = 1:length(ctypes)
            delx = diff(ceps{ir});
            if p.expandStart && delx <=5
                ceps{ir} = [ceps{ir}(1), ceps{ir}(1) + p.expandStartMagn*delx];
            end
            myRectangle(cHand(zi), ceps{ir},...
                [cy(1), cy(1) + cheight], ctypes(ir), [0 0 0],[1 1 1], [.5 .5 .5])
        end
    end
end
end

function out = plotTimeCourse(data, p, cax)
cols = p.cols;
d    = p.d;
n    = p.n;
fs   = p.fs;

%% B: plot IS data
% Add plots
% cla(sp.B);
% sp.B.PositionConstraint = 'innerposition';
% ylimits.fin  = [-.85 .85];
delGR  = [-5 5];
for ig = 1:n.groups % Plot in reverse order to overlap curves conveniently
    %     scatter(1:n.strTot, summary.m{ig}, 20,  gr.cols{ig}, 'filled');
    hold on
    if strcmp(p.type,'data')
        [lh(ig), shadh(ig)] = boundedline(1:p.n.totStrides, data.m{ig}, data.se{ig}, ...
            '-o', 'nan', 'gap');
        shadh(ig).FaceColor    = cols.groups(ig,:);
        shadh(ig).FaceAlpha = .5;
    elseif strcmp(p.type,'model')
        lh(ig) = plot(nanTrans(-p.parad{ig}.y,p.parad{ig}.transitions),'color',...
            cols.groups(ig,:));
        if ig==2 
            indBS = 1:150;
            plot(indBS, -p.parad{1}.y(indBS),'color',...
            cols.groups(1,:),'LineStyle','--','LineWidth', 1);
        
            indA1 = 152:750;
            plot(indA1, -p.parad{1}.y(indA1),'color',...
            cols.groups(1,:),'LineStyle','--','LineWidth', 1);
        end
        indAfTransitions    = p.parad{ig}.transitions + 1;
        iat    = indAfTransitions;
        iat(1) = iat(1) + delGR(ig);
        scatter(iat, -p.parad{ig}.y(indAfTransitions),p.dotSize,cols.groups(ig,:),'filled');
    end
    lh(ig).LineWidth = 1;
    lh(ig).Color           = cols.groups(ig,:);
    lh(ig).MarkerFaceColor = cols.groups(ig,:);
    axis tight
end
% ylimits.init = ylim();
ylim(p.ylimits);

% Add patches
p.patchParams.sav.ylims = p.ylimits;
p.patchParams.int.ylims = p.ylimits;
addParadPatches(cax, p.parad{1}, p.patchParams.sav);
addParadPatches(cax, p.parad{1}, p.patchParams.int);

% Add epoch names & rectangles
pr   = struct('expandStart', p.expandStart, 'expandStartMagn', p.expandStartMagn, 'height', p.rect);
rect = addEpochsNames(p.epochs, -d, p.ylimits, 10, cols.epochs, pr, p.addEpochLabels);
% hiddenAx = myHiddenAxes(cax.Position, handles.top);

if strcmp(p.type,'model')
    ylab(2) = ylabel('Step Length Asymm.', 'FontSize', fs.xxs);
    ylab(2).Position = ylab(2).Position + [0 -d*1.1 0] ;
end
% if strcmp(p.type,'data')
%     ylab(2) = ylabel('Step Length Asymm.', 'FontSize', fs.xxs);
%     ylab(2).Position = ylab(2).Position + 0*[0 -d 0] ;
% end

% xlabel('Strides','FontSize', fs.xxs)
% text(mean(cax.XLim), ylimits.fin(1)*1.4, 'Strides', 'FontSize', fs.xxs, 'HorizontalAlignment', 'center');
xticks(p.parad{1}.transitions);
cax.XTickLabel{4} = ['      ' cax.XTickLabel{4}]; %Shift to the right
out.rect = rect;
end

function rect = addEpochsNames(ep, dataLim, axesLim, fsize, epcols, pr, addEpochLabels)
cax = gca;
ybottom = axesLim(1);
neps    = length(ep.def);
for e = 1:neps
    ind.Start   = [1 4 6 9];
    ind.Early   = [2 7];
    ind.End     = [3 5 8 10];
    %     masks.Start = mySetTrue(ind.Start,neps);
    %     masks.Early = mySetTrue(ind.Start,neps);
    %     masks.End   = mySetTrue(ind.Start,neps);
    
    ccol = epcols(e,:);
    delLine = .005; %Delta epoch's line
    
    
    x = ep.def{e};
    
    
    % Add rectangle
    yRec.top   = dataLim(1);
    if isempty(pr.height)
        yRec.bot   = dataLim(1) - 13*delLine;
    else
        rHeight    = pr.height.fNormHeight*diff(axesLim)/cax.Position(4);
        yRec.bot   = dataLim(1) - rHeight;
    end
    xRec.left  = x(1);
    xRec.right = x(2);
    %     dx = 3; % To make rectangles larger
    switch e
        case num2cell(ind.Start)
            type = 1;
            if pr.expandStart
                xRec.right = xRec.right + pr.expandStartMagn*(xRec.right - xRec.left);
            end
        case num2cell(ind.Early)
            type = 2;
            %             xRec.right = xRec.right + dx*(xRec.right - xRec.left);
        case num2cell(ind.End)
            type = 3;
            %             xRec.left = xRec.left - dx*(xRec.right - xRec.left);
    end
    myRectangle(cax, xRec, yRec, type, [0 0 0], [1 1 1], [.5 .5 .5]);
    
    % Previous end
    %     if e>1 && e <= 3
    %         pExtent = th(e-1).Extent;
    %         pXend   = sum(pExtent([1 3]));
    %     end
    switch e
        case 1
            xpos = x(2);
            horAlignment = 'right';
        case 2
            xpos = x(1);
            horAlignment = 'left';
        case 3
            xpos = x(2);
            horAlignment = 'right';
    end
    
    ystring = ybottom;
    vertAli = 'baseline';
    
    if addEpochLabels
        switch e
            case num2cell([1 2 3])
                hold on, th(e) = text(xpos, ystring, ep.names{e},'FontSize', fsize,...
                    'color', [0 0 0], 'VerticalAlignment', vertAli , 'HorizontalAlignment', horAlignment);
        end
    end
    % Add vertical line from epoch's definition to its string
    %     if ~any(e==el1) %Row 2
    %         hold on, line(mean(x)*ones(2,1), [dataLim(1) - delLine, dataLim(1) - delyString], 'Color', ccol, 'LineWidth', 1)
    %     end
    %       a = annotation('textarrow',x,y,'String','y = x ');
end
% Compute normalized height of rectangle
% myRectangle(gca, xRec, yRec, type, [0 0 0], [1 1 1], [.5 .5 .5]);
rect.height      = yRec.top-yRec.bot; % Height
rect.range       = diff(cax.YLim);
rect.normHeight  = cax.Position(4)*rect.height/rect.range; %Height normalized with respect to figure
end

function  myRectangle(cax, xr, yr, type, colPat, colStr, colEdge)
if isstruct(xr)
    x = [xr.left xr.right];
else
    x = xr;
end
if isstruct(yr)
    y = [yr.bot  yr.top];
else
    y = yr;
end
% yhline.top = dataLim(1);
% yhline.bot = dataLim(1) - 13*delLine;
%     yhline.mid = dataLim(1) - 10*delLine;
%     hold on, line(x, yhline*ones(1,2), 'Color',ccol,'LineStyle','-', 'LineWidth', 4)
hold on
patch(cax, [x fliplr(x)],[y(1)*ones(1,2) y(2)*ones(1,2)], colPat)
hold on

switch type
    case 1
        line(cax, x, y,'Color', colStr,'linewidth',1);
    case 2
        line(cax, x, fliplr(y),'Color', colStr,'linewidth',1);
    case 3
        line(cax, x, y,'Color', [1 1 1],'linewidth',1);
        line(cax, x, fliplr(y),'Color', colStr,'linewidth',1);
end
patch(cax, [x fliplr(x)],[y(1)*ones(1,2) y(2)*ones(1,2)],...
    colPat, 'Facecolor', 'none', 'Linewidth', .5, ...
    'EdgeColor', colEdge);
end