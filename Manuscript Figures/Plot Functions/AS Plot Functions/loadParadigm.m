function [paradigm, epOrder] = loadParadigm(paradName, varargin)
if nargin>1
    d = varargin{1};
else
    d = 1;
end
switch paradName
    case 'ANA'
        n            = struct('epochs',4,'trialsPerEpoch',[150,600,600,600]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','NL','A2'},...
                              'printID', {'N','A_1','N_L','A_2'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = d*ones(1,n.trialsPerEpoch(4));
    case 'SHAB'
        n            = struct('epochs',4,'trialsPerEpoch',[150,120,600,600]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','NL','A2'},...
                              'printID', {'N','A_1','N_L','A_2'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = d*ones(1,n.trialsPerEpoch(4));
    case 'GR'
        n            = struct('epochs',4,'trialsPerEpoch',[150,600,600,600]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','NL','A2'},...
            'printID', {'N','A_1','N_L','A_2'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = [zeros(1,tclamp),linspace(0,d,n.trialsPerEpoch(2)-2*tclamp),zeros(1,tclamp)];
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = d*ones(1,n.trialsPerEpoch(4));
    case 'SAV'
        n            = struct('epochs',6,'trialsPerEpoch',[150,600,1200,150,600,150]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','NL','NS','A2','WO'},...
                              'printID', {'N','A_1','N_L','N_S','A_2','W_O'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = zeros(1,n.trialsPerEpoch(4));
        epochs(5).u  = d*ones(1,n.trialsPerEpoch(5));
        epochs(6).u  = zeros(1,n.trialsPerEpoch(6));
    case 'INT'
        n            = struct('epochs',6,'trialsPerEpoch',[150,600,1200,150,600,150]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','B','NS','A2','WO'},...
            'printID', {'N','A_1','B','N_S','A_2','W_O'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = [-d*ones(1,.5*n.trialsPerEpoch(3)) linspace(-d,0,.5*n.trialsPerEpoch(3))];
        epochs(4).u  = zeros(1,n.trialsPerEpoch(4));
        epochs(5).u  = d*ones(1,n.trialsPerEpoch(5));
        epochs(6).u  = zeros(1,n.trialsPerEpoch(6));
    otherwise
        warning('unrecognized paradigm')
end
        n.ttot = sum(n.trialsPerEpoch);
        t      = 1:n.ttot;
        % Store epoch start
        e = 1;
        epochs(e).tstart = 1;
        epochs(e).tend   = epochs(e).tstart + n.trialsPerEpoch(e) - 1;
        midEpochs(e)     = epochs(e).tstart + round(.5*n.trialsPerEpoch(e) - 1);
        for e = 2:n.epochs
            epochs(e).tstart = epochs(e-1).tstart + n.trialsPerEpoch(e-1);
            epochs(e).tend   = epochs(e).tstart   + n.trialsPerEpoch(e) - 1;
            transitions(e-1) = epochs(e).tstart; 
            midEpochs(e)     = epochs(e).tstart + round(.5*n.trialsPerEpoch(e) - 1); 
        end
        
        % Concatenate u
        u = [epochs(:).u];
        epOrder  = cell2struct(num2cell(1:n.epochs)',{epochs.id}');
        paradigm = struct('epochs',epochs,'transitions',transitions,...
                          'midEpochs',midEpochs,'n',n,...
                          'u',u,'t',t);
end