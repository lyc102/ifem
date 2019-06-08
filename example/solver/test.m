[node,elem] = squaremesh([0 1 0 1],0.5);
for k = 1:6
  [node,elem] = uniformrefine(node,elem);
end
pde.f = 1;
pde.g_D = 0;
option.solver = 'none';
[soln,eqn] = PoissonP3(node,elem,[],pde,option);
fprintf('Number of unknowns: %8.0u\n',length(eqn.b))
tic; disp('Direct solver'); x1 = eqn.A\eqn.b; toc;
tic; x2 = mg(eqn.A,eqn.b,elem); toc;
tic; x3 = amg(eqn.A,eqn.b); toc;
format shorte
fprintf('Difference between direct and mg, amg solvers %0.2g, %0.2g \n',...
         norm(x1-x2)/norm(eqn.b),norm(x1-x3)/norm(eqn.b));