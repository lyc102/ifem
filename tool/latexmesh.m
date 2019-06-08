function latexmesh(node,elem,filename,scaling,landscape)
%% LATEXMESH draw a mesh using pstrick
%  run 'latex mesh.tex' to get the figure
%
% load lakemesh.mat
% latexmesh(node,elem,'Lshape.tex');
% latexmesh(node,elem);
% 
% If the file name is not given, the default one is 'mesh.tex'.
%
% Copyright Long Chen.

if ~exist('scaling','var'), scaling = 2; end
if ~exist('landscape','var'), landscape = 'landscape'; end

%% Rescale the corrdinate
maxxy = max(node(node(:,1)~=Inf,:));
minxy = min(node(node(:,1)~=Inf,:));
node(:,1) = node(:,1)-0.5*(minxy(1)+maxxy(1));
node(:,2) = node(:,2)-0.5*(minxy(2)+maxxy(2));
node = scaling*node/(0.5*min(maxxy-minxy));
shift = min(node(node(:,1)~=Inf,:));
node(:,1) =  node(:,1)  - 5/16*shift(:,1);
node(:,2) =  node(:,2)  - 5/16*shift(:,2);
box = ceil(max(node(node(:,1)~=Inf,:)));

%% Draw the mesh using pstrick
if (nargin <=2) || isempty(filename), filename='mesh.tex'; end
fid = fopen(filename,'wt');
if (strcmp(landscape,'landscape') == true)
    fprintf(fid,'\\documentclass[landscape,10pt]{article} \n');
else
    fprintf(fid,'\\documentclass[10pt]{article} \n');
end
fprintf(fid,'\\usepackage{times,amsmath,amsbsy,amssymb} \n');
fprintf(fid,'\\usepackage{graphicx,epic,pstricks} \n');
fprintf(fid,'\\begin{document} \n');
fprintf(fid,'\\pagenumbering{gobble} \n');
fprintf(fid,'\\begin{figure} \n');
fprintf(fid,'\\begin{center} \n');
fprintf(fid,'\\psset{unit=1cm,linewidth=0.5pt} \n');
fprintf(fid,['\\begin{pspicture}(' int2str(box(1)) ',' int2str(box(2)) ')\n']);
% output grid
for t=1:size(elem,1)
    if size(elem,2) == 3 % triangular grids
    fprintf(fid,'\\pspolygon(%2.3f,%2.3f)(%2.3f,%2.3f)(%2.3f,%2.3f) \n',...
       node(elem(t,1),1),node(elem(t,1),2),node(elem(t,2),1),node(elem(t,2),2),...
       node(elem(t,3),1),node(elem(t,3),2));
    elseif size(elem,2) == 4 % rectangular grids
    fprintf(fid,'\\pspolygon(%2.3f,%2.3f)(%2.3f,%2.3f)(%2.3f,%2.3f)(%2.3f,%2.3f) \n',...
       node(elem(t,1),1),node(elem(t,1),2),node(elem(t,2),1),node(elem(t,2),2),...
       node(elem(t,3),1),node(elem(t,3),2),node(elem(t,4),1),node(elem(t,4),2));
    end        
end
fprintf(fid,'\\end{pspicture} \n');
% fprintf(fid,'\\caption{Mesh} \n');
fprintf(fid,'\\end{center} \n');
fprintf(fid,'\\end{figure} \n');
fprintf(fid,'\\end{document}  \n');
fclose(fid);