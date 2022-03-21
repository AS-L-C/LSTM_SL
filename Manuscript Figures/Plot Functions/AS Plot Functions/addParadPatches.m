function addParadPatches(sph, parad, params)
% IDs      = {parad.epochs.printID};
xlefts     = [parad.epochs.tstart] - params.dx;
xrights    = [parad.epochs.tend]   + params.dx;

% ytops    = params.ylims(2)*ones(1,ncs);
% ybots    = params.ylims(1)*ones(1,ncs);
patchInds  = params.patchInd;
patchPos   = params.patchPos;
nPatches   = length(patchInds);
assert(nPatches==length(patchPos));
cols       = params.colors;
xlefts     = xlefts(patchInds);
xrights    = xrights(patchInds);
for ic = 1:nPatches % For every condition
        switch patchPos(ic)
            case -1
                ytop = 0;
                ybot = params.ylims(1);
            case  0
                ytop = params.ylims(2);
                ybot = params.ylims(1);
            case  1
                ytop = params.ylims(2);
                ybot = 0;
        end
        addPatches(sph, cols(ic,:), xlefts(ic), xrights(ic), ybot, ytop,...
                                                     0, 0, params.alpha); 
end



end