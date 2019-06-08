%% Check Transfer Operator for edge element in Three dimensions
%
% This script is to check transfer operators, prolongation and restriction,
% for the edge element in three dimensions.
%
% Written by Jie Zhou.

%% Test data
%u  = inline('[1-x(:,2)-x(:,3), x(:,1), x(:,1)]','x');
%u  = inline('[x(:,2), 1-x(:,1)-x(:,3), x(:,2)]','x');
u   = inline('[x(:,3), x(:,3), 1-x(:,1)-x(:,2)]','x');
%u  = inline('[-x(:,2), x(:,1), 0*x(:,1)]','x');
%u = inline('[-x(:,3), 0*x(:,1), x(:,1)]','x');
%u  = inline('[0*x(:,3), -x(:,3), x(:,2)]','x');
%u = inline('[x(:,3)+4, x(:,3), 9-x(:,1)-x(:,2)]','x');
%% Test 1: Uniform refine
[node,elem] = cubemesh([0,1,0,1,0,1],1.0/16);
[nodef,elemf] = uniformrefine3(node,elem);
[~,edgef] = dof3edge(elemf);
[~,edge]  = dof3edge(elem);
pro   = transferedgered3(elem,elemf);
%Test two different uI by interpolate and  prolongation
uI_c  = edgeinterpolate(u,node,edge);
uI_f  = edgeinterpolate(u,nodef,edgef);
u_c2f = pro*uI_c;
%pro
disp('Check the edgetransfer operator is correct or not:')
disp(norm(uI_f-u_c2f))

%% Test 2: Uniform coarsen
[nodef,elemf] = cubemesh([0,1,0,1,0,1],1.0/16);
elem = uniformcoarsen3red(elemf);
[~,edgef] = dof3edge(elemf);
[~,edge]  = dof3edge(elem);
pro   = transferedgered3(elem,elemf);
%Test two different uI by interpolate and  prolongation
uI_c  = edgeinterpolate(u,node,edge);
uI_f  = edgeinterpolate(u,nodef,edgef);
u_c2f = pro*uI_c;
%pro
disp('Check the edgetransfer operator is correct or not:')
disp(norm(uI_f-u_c2f))