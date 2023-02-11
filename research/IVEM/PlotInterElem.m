function PlotInterElem(vert,intpt)


scatter3(vert(1,1),vert(1,2),vert(1,3),'k');
hold on

scatter3(vert(2,1),vert(2,2),vert(2,3),'k');
scatter3(vert(3,1),vert(3,2),vert(3,3),'k');
scatter3(vert(4,1),vert(4,2),vert(4,3),'k');

if size(intpt,1) == 3
    scatter3(intpt(1,1),intpt(1,2),intpt(1,3),'r');
    scatter3(intpt(2,1),intpt(2,2),intpt(2,3),'b');
    scatter3(intpt(3,1),intpt(3,2),intpt(3,3),'m');
elseif size(intpt,1) == 4
    scatter3(intpt(1,1),intpt(1,2),intpt(1,3),'k');
    scatter3(intpt(2,1),intpt(2,2),intpt(2,3),'b');
    scatter3(intpt(3,1),intpt(3,2),intpt(3,3),'k');
    scatter3(intpt(4,1),intpt(4,2),intpt(4,3),'y');
end

axis equal

hold off
