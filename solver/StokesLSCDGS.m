function [u,p] = StokesLSCDGS(u,p,f,g,A,B,auxMat,elem,smootherOpt,Ai,Si,SSi,Res,Pro)
%% STOKESLSCDGSTRI Matrix-vecotr form DGS relaxation for Stokes eqns
%
%************************** Algorithm Description *************************
% 
%  In matrix form, we need to solve the equations
%
%         Lx =  |A B' | |u| = |f|
%               |B -D | |p| = |0|
%  Distributive matrix will be given by
%
%         M = |I            B'|
%             |0 -inv(BB')BAB'|,
%
%  Therefore, the transformed matrix T will be
%
%         T = L*M = |A  PAB'               |    with P = I-B'inv(BB')B
%                   |B  BB'+D*inv(BB')BAB' |
%
%         tildeT = |Su   0|    with Su and Sp smoother for A and BB'.
%                  |B   Sp|
%
%  The matrix form DGS update can be written as
%
%    |uk+1|   |uk|                 |ru|     with ru = f-Auk -B'pk,
%    |    | = |  | + M*inv(tildeT)*|  |
%    |pk+1|   |pk|                 |rp|     and  rp = 0-Buk
%
%********************** End-of- Algorithm Description *********************
%
%  Created by Ming Wang (with discussion with Long Chen) at Jan, 2012.

%% Set up of smoothers
if isfield(smootherOpt,'smoothingstep')
    smoothingstep = smootherOpt.smoothingstep;
else
    smoothingstep = 2;
end
if isfield(smootherOpt,'smootherSp')
    smootherSp = upper(smootherOpt.smootherSp);
else
    smootherSp = 'SGS';
end
if isfield(smootherOpt,'smootherbarSp')
    smootherbarSp = upper(smootherOpt.smootherbarSp);
else
    smootherbarSp = 'SGS';
end
if isfield(smootherOpt,'smootherbarSpPara')
    smootherbarSpPara = smootherOpt.smootherbarSpPara;
else
    smootherbarSpPara = 1;
end
if strcmp(smootherbarSp,'VCYCLE') || strcmp(smootherSp,'VCYCLE')
    % mg Parameters for barSp
    optionmg.solvermaxit = 1;
%     optionmg.tol = 0.1;
    optionmg.solver = 'VCYCLE'; 
    optionmg.smoothingstep = 2;
    optionmg.printlevel=0; 
    optionmg.setupflag = 0;
end

%% initialize smoother 
Bt = auxMat.Bt;
BBt = auxMat.BBt;
BABt = auxMat.BABt;
Su = auxMat.Su;
Sp = auxMat.Sp;
Spt = auxMat.Spt;
DSp = auxMat.DSp; 

%% DGS relaxation step
for k = 1: smoothingstep
    % Step 1: relax Momentum eqns
    u = u + Su\(f-Bt*p-A*u);
    % Step 2: relax transformed Continuity eqns
    rp = g - B*u;
    switch(smootherSp)
        case 'SGS'
            dq = Spt\(DSp.*(Sp\rp)); % symmetric Gauss-Seidel iteration
        case 'GS'
            dq = Sp\rp; % Gauss-Seidel iteration
        case 'VCYCLE'
            if exist('elem','var')
                dq = mg(BBt,rp,elem,optionmg,Ai,Si,SSi,Res,Pro); % mg vcycle
            else
                dq = amg(BBt,rp,optionmg); % mg vcycle
            end
    end
    % Step 3: transform the correction back to the original variables
    %   Step 3.1: velocity
    u = u + Bt*dq;
    %   Step 3.2: pressure
    dq = BABt*dq;
%     dq = dq - mean(dq);
    switch(smootherbarSp)
        case 'GS'
            dp = Spt\dq;
        case 'SGS'
            dp = Sp\(DSp.*(Spt\dq));
%             dp = dp + Sp\(DSp.*(Spt\(dq - BBt*dp)));
        case 'VCYCLE'
%            dp = BBt\dq;
            if exist('elem','var')
%                 dq =dq-mean(dq);
                dp = mg(BBt,dq,elem,optionmg,Ai,Si,SSi,Res,Pro); % mg vcycle;
            else
                dp = amg(BBt,dq,optionmg);
            end
    end
    p = p - smootherbarSpPara*dp;     
end