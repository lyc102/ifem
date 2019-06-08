function option = femoption(option)
%% FEMOPTION provides default options for femPoisson
%
% option = femoption(option)
%
%  - options.elemType:  type of finite element
%  - options.maxIt      
%  - options.maxN
%  - options.L0
%  - options.refType    mesh refinement type 
%  - options.printlevel
%  - options.plotflag
%  - options.rateflag
%  - options.dispflag
%  - options.tol
%  - options.lumpflag


if ~isfield(option,'elemType')
    option.elemType = 'P1';    
end

if ~isfield(option,'maxIt')
    option.maxIt = 4;    
end

if ~isfield(option,'maxN')
    option.maxN = 2e5;    
end

if ~isfield(option,'L0')
    option.L0 = 0;    
end

if ~isfield(option,'refType')
    option.refType = 'red';
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

if ~isfield(option,'tol')
    option.tol = 1e-8;    
end

if ~isfield(option,'lumpflag') % mass lumping
    option.lumpflag = 0;       % no mass lumping
end