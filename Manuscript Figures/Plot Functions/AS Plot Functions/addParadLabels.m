function addParadLabels(sph, parad, fsize, ylims, posStrings, indStrings)

nstrings = length(indStrings);
assert(nstrings==length(posStrings));
printIDs   = {parad.epochs.printID};
xlefts     = [parad.epochs.tstart]; 
xrights    = [parad.epochs.tend];   
xs         = round(.5*xlefts+.5*xrights);
% ytops    = params.ylims(2)*ones(1,ncs);
% ybots    = params.ylims(1)*ones(1,ncs);

for is = 1:nstrings
    switch posStrings(is)
        case -1
            y = mean([ylims(1) 0]);
        case 0 
            y = 0;
        case 1
            y = mean([0 ylims(2)]);
    end
    text(sph, xs(indStrings(is)), y, printIDs{indStrings(is)},...
           'FontSize', fsize, 'HorizontalAlignment', 'center');
end


end