function option = femStokesoption(option)

if ~isfield(option,'elemType')
    option.elemType = 'P2P1';
end

if strcmp(option.elemType,'P1') % default for Poisson
    option.elemType = 'P2P1';   % modify for Stokes
end

if ~isfield(option,'smoothingstep')
    option.smoothingstep = 2;
end

if ~isfield(option,'smootherSp')
    option.smootherSp = 'SGS';
end

if ~isfield(option,'smootherbarSp')
    option.smootherbarSp = 'SGS';
end

if ~isfield(option,'smootherType')
    option.smootherType = 'LSCDGS';
end