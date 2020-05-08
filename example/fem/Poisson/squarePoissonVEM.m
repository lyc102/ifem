clear all
clc
close all
%% Load the information of mesh and function 
dirName = 'data/';
%load the function
 pde = vemcosdata;
%pde = vemdata;
%load the elements
%nameV = [128,512,2048,8192,32768];
nameV = [16,64,256,1024,4096,20014];
%nameV = 131072;
errorH1 = zeros(1,length(nameV)); %initialize the error
errorMax = zeros(1,length(nameV)); %initialize the error
solverTime = zeros(length(nameV),1); 
assembleTime = zeros(length(nameV),1);
N = zeros(1,length(nameV)); %initialize the error
for k = 1:length(nameV)
    % load mesh
    aa = [dirName,'E',num2str(nameV(k)),'.mat'];
    load(aa);
    
    %% The exact value and right hand side function for Possion equation
    % get the approximate value u and corresponding error
    tic;
    [u,A,asst,solt] = PoissonVEM(Node,Element,pde);
    assembleTime(k) = asst;
    solverTime(k) = solt;
    errorH1(k) = sqrt((u-pde.u(Node))'*A*(u-pde.u(Node))); % get the error in H1 norm
    errorMax(k) = max(abs(u-pde.u(Node)));
    N(k) = length(u);
end
N = N';
nameV = nameV';
errorH1 = errorH1';
errorMax = errorMax';

%% Plot convergence rates
figure;
set(gcf,'Units','normal');
set(gcf,'Position',[0.25,0.25,0.55,0.4]);
showrate2(N(1:k),errorH1(1:k),3,'r-+','||u_I-u_h||_A',...
    N(1:k),errorMax(1:k),3,'b-+','||u_I-u_h||_{\infty}');


%% Output
err = struct('N',N,'elem',nameV,'H1',errorH1(1:k),'uIuhMax',errorMax(1:k));
time = struct('N',N,'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k));
% %% Display error and time
option.dispflag =1;
 if option.dispflag
     display('Table: Error')
   colname = {'#Dof','# of elements','||DuI-Du_h||','||uI-u_h||_{max}'};     
   %disptable(colname,err.N,[],err.h,'%0.2e',err.L2,'%0.5e',err.H1,'%0.5e',err.uIuhH1,'%0.5e');
     disptable(colname,err.N,[],err.elem,[],err.H1,'%0.5e',err.uIuhMax,'%0.5e');
     display('Table: CPU time')
    colname = {'#Dof','Assemble','Solve'};
    disptable(colname,time.N,[],time.assemble,[],time.solver,[]);
%     display('Table: CPU time')
%      colname = {'#Dof','Assemble','Solve','Error','Mesh'};
%     disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
%         time.err,'%0.2e',time.mesh,'%0.2e');
 end
