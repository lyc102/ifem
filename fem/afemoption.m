function option = afemoption(option,dim)

if ~isfield(option,'viewcut')
    option.viewcut = [];
end

if ~isfield(option,'viewangle')
    switch dim
        case 2
            option.viewangle = [-50,12];
        case 3
            option.viewangle = [74 20];
    end
end

if ~isfield(option,'elemType')
    option.elemType = 'P1';    
end

if ~isfield(option,'maxIt')
    option.maxIt = 50;    
end

if ~isfield(option,'maxN')
    switch dim
        case 2
            option.maxN = 8e3;    
        case 3
            option.maxN = 1e4;
    end
end

if ~isfield(option,'theta')
    option.theta = 0.35;    
end

if ~isfield(option,'L0')
    option.L0 = 0;    
end

if ~isfield(option,'printlevel')
    option.printlevel = 1;    
end

if ~isfield(option,'plotflag')
    option.plotflag = 1;    
end

if ~isfield(option,'rateflag')
    option.rateflag = 1;    
end

if ~isfield(option,'dispflag')
    option.dispflag = 1;    
end

if ~isfield(option,'refType')
    option.refType = 'bisect';    
end

if ~isfield(option,'estType')
    option.estType = 'residual';    
end

if ~isfield(option,'markType')
    option.markType = 'L2';    
end

if ~isfield(option,'coarsenflag')
    option.coarsenflag = 0;    
end