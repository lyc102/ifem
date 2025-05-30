classdef testtet2hex < matlab.unittest.TestCase
    % Test class for tet2hex conversion functionality
    
    methods(Test)
        function testUniformGrid(testCase)
            % Test tet2hex conversion on uniform grid
            
            % Create uniform grid
            [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
            
            % Convert tetrahedra to hexahedra
            [node_hex,hexelem] = tet2hex(node,elem,HB);
            
            % Verify output dimensions
            testCase.verifyEqual(size(node_hex,2), 3, ...
                'Node coordinates should be 3D');
            testCase.verifyGreaterThan(size(hexelem,1), 0, ...
                'Should generate at least one hexahedral element');
            testCase.verifyEqual(size(hexelem,2), 8, ...
                'Hexahedral elements should have 8 nodes');
            
            % Verify all node indices are valid
            testCase.verifyLessThanOrEqual(max(hexelem(:)), size(node_hex,1), ...
                'All node indices should be within valid range');
            testCase.verifyGreaterThanOrEqual(min(hexelem(:)), 1, ...
                'All node indices should be positive');
        end
        
        function testAdaptiveGrid(testCase)
            % Test tet2hex conversion on adaptively refined grid
            
            % Create initial grid
            [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1/4);
            [node,elem] = delmesh(node,elem,'x<0 & y<0');
            bdFlag = setboundary3(node,elem,'Neumann');
            
            % Perform adaptive refinement
            for k = 1:2
                eta = abs(sign(testCase.sphereFunction(node(elem(:,1),:))) + ...
                         sign(testCase.sphereFunction(node(elem(:,2),:))) + ...
                         sign(testCase.sphereFunction(node(elem(:,3),:))) + ...
                         sign(testCase.sphereFunction(node(elem(:,4),:))));
                refineElem = find(eta < 4);
                [node,elem,bdFlag,HB] = bisect3(node,elem,refineElem,bdFlag,HB);
            end
            
            % Convert to hexahedral mesh
            [node_hex,hexelem] = tet2hex(node,elem,HB);
            
            % Verify output dimensions
            testCase.verifyEqual(size(node_hex,2), 3, ...
                'Node coordinates should be 3D');
            testCase.verifyGreaterThan(size(hexelem,1), 0, ...
                'Should generate at least one hexahedral element');
            testCase.verifyEqual(size(hexelem,2), 8, ...
                'Hexahedral elements should have 8 nodes');
            
            % Verify mesh quality
            testCase.verifyLessThanOrEqual(max(hexelem(:)), size(node_hex,1), ...
                'All node indices should be within valid range');
            testCase.verifyGreaterThanOrEqual(min(hexelem(:)), 1, ...
                'All node indices should be positive');
        end
        
        function testBoundaryHandling(testCase)
            % Test that boundary conditions are handled properly
            
            [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1/2);
            bdFlag = setboundary3(node,elem,'Neumann');
            
            % Find boundary faces before conversion
            [~,bdFace] = findboundary3(elem,bdFlag);
            testCase.verifyGreaterThan(size(bdFace,1), 0, ...
                'Should have boundary faces');
            
            % Convert to hex mesh
            [node_hex,hexelem] = tet2hex(node,elem,HB);
            
            % Verify conversion doesn't break
            testCase.verifyEqual(size(node_hex,2), 3, ...
                'Conversion should preserve 3D coordinates');
        end
        
        function testVisualization(testCase)
            % Test that visualization functions work without errors
            
            [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
            
            % Test mesh visualization
            testCase.verifyWarningFree(@() showmesh3(node,elem), ...
                'showmesh3 should execute without warnings');
            
            % Test node finding
            testCase.verifyWarningFree(@() findnode3(node,'all'), ...
                'findnode3 should execute without warnings');
        end
    end
    
    methods(Static)
        function s = sphereFunction(p)
            % Sphere level set function for adaptive refinement
            s = sum(p.^2,2) - (0.5)^2;
        end
    end
end