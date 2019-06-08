function plot_error_table(dof1,error1,dof2,error2,dof3,error3,dof4,error4)

set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',1,...
      'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
      'defaulttextfontsize',18);

mindof = min([min(dof1) min(dof2) min(dof3) min(dof4)]);
maxdof = max([max(dof1) max(dof2) max(dof3) max(dof4)]);
% minerr = min([min(error1) min(error2) min(error3) min(error4)]);
maxerr = max([max(error1) max(error2) max(error3) max(error4)]);

r = -1/3;  
s = 0.42*maxerr/mindof^r;  
loglog(dof1,s*dof1.^r,'r','LineWidth',2.5);
hold on;
loglog(dof1,error1,'-.c+','MarkerSize',13,'LineWidth',2.5);
hold on;
loglog(dof2,error2,'-.ob','MarkerSize',13,'LineWidth',2.5);
hold on;
loglog(dof3,error3,'-.kv','MarkerSize',13,'LineWidth',2.5)
hold on
loglog(dof4,error4,'-.md','MarkerSize',13,'LineWidth',2.5) 

grid on
% figure (1);
xlabel('Degrees of Freedom (DOFs)','interpreter','latex', 'FontSize', 22)
ylabel('Error','interpreter','latex', 'FontSize', 22)
               
h = legend('$$\mathrm{ DOFs}^{\mathrm{-}\frac{1}{3}}$$', ...
'$$s = 0.2$$', ...
'$$s = 0.4$$', ...
'$$s = 0.6$$', ...
'$$s = 0.8$$', ...
'Location','NorthEast' );
set(h,'Interpreter','latex','fontsize',22);

axis([0.75*mindof 1.25*maxdof 0.8*s*maxdof^r 1.2*maxerr])

%since the Latex labels are invisble for plotting
%and printing purposes, we have to create a margin for the plot
%so that we can actually read the labels
set(gcf,'Units','normal')
set(gca,'Position',[0.14 0.17 0.82 0.75])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.0 0.0 7 5]);
%Set the paper to have width 5 and height 5.'$$\textrm{Estimator(s=0.2)}$$', ...
set(gcf, 'PaperSize', [7 5]); 