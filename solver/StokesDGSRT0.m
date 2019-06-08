function [u,p] = StokesDGSRT0(A,B,u,p,f,g,itStep,auxMat)
%% STOKESDGS
%  Created by Ming Wang. Discussed with Long Chen. 
% Reference: 

Su = auxMat.trilA; 
Bt = auxMat.Bt;
Sp = auxMat.Sp;
invMp = auxMat.invMp;

%% DGS relaxation step
for k = 1: itStep
    % Step 1: relax Momentum eqns twice
    for i = 1:2
        u = u + Su\(f-Bt*p-A*u);
    end
    % Step 2: relax Continuity eqns
    rp = g - B*u;
%     rp = rp - mean(rp);  Whether need this step ?????
    dq = Sp\rp;
    % Step 3: update u and p
    u = u + Bt*dq;
    p = p - invMp*rp;
end