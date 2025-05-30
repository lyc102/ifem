classdef testnodeinterpolate < matlab.unittest.TestCase
     % Test class for node interpolation functionality
     
     methods(Test)
         function testNodeInterpolationWithBisection(testCase)
             % Test node interpolation through multiple bisection operations
             
             % Setup initial mesh
             [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],2);
             u = ones(size(node,1),1);
             
             % Perform multiple bisection operations
             [node,elem,~,HB] = bisect3(node,elem,[1 2],[],HB);
             [node,elem,~,HB] = bisect3(node,elem,'all',[],HB);
             [node,elem,~,HB] = bisect3(node,elem,[1 2 3],[],HB);
             [node,elem,~,HB] = bisect3(node,elem,[1 2 3],[],HB);
             
             % Test interpolation after refinement
             u_refined = nodeinterpolate(u,HB);
             
             % Verify interpolation preserves constant function
             testCase.verifyEqual(u_refined, ones(size(node,1),1), ...
                 'AbsTol', 1e-14, ...
                 'Interpolation should preserve constant function');
             
             % Test coarsening and interpolation
             [node_coarse,elem_coarse,~,HB_coarse,indexMap] = coarsen3(node,elem,'all',[],HB);
             u_coarse = nodeinterpolate(u_refined,indexMap);
             
             % Verify dimensions match after coarsening
             testCase.verifyEqual(size(u_coarse,1), size(node_coarse,1), ...
                 'Solution size should match coarsened mesh size');
             
             % Verify interpolation preserves constant function after coarsening
             testCase.verifyEqual(u_coarse, ones(size(node_coarse,1),1), ...
                 'AbsTol', 1e-14, ...
                 'Interpolation should preserve constant function after coarsening');
         end
         
         function testMeshVisualization(testCase)
             % Test that mesh visualization doesn't throw errors
             [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],2);
             [node,elem,~,HB] = bisect3(node,elem,[1 2],[],HB);
             
             % Test that showmesh3 can be called without error
             testCase.verifyWarningFree(@() showmesh3(node,elem), ...
                 'showmesh3 should execute without warnings');
         end
     end
 end