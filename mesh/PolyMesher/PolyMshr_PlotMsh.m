function PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load)
clf; axis equal; axis off; hold on;
Element = Element(1:NElem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
if exist('Supp','var')&&~isempty(Supp) %Plot Supp BC if specified
  plot(Node(Supp(:,1),1),Node(Supp(:,1),2),'b>','MarkerSize',8);
end
if exist('Load','var')&&~isempty(Load) %Plot Load BC if specified
  plot(Node(Load(:,1),1),Node(Load(:,1),2),'m^','MarkerSize',8);
end