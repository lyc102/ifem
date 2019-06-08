% AFEM
%
% Files
%   bisect            - bisect a 2-D triangulation.
%   bisect3           - bisect a 3-D triangulation.
%   coarsen           - coarsen a 2-D triangulation.
%   coarsen3          - coarse a 3-D triangulation.
%   edgeinterpolate   - interpolate to the lowest order edge finite element space.
%   edgeinterpolate1  - interpolate to the linear edge finite element space.
%   edgeinterpolate2  - interpolate to the quadratic (1st type) edge finite element space.
%   eleminterpolate   - interpolate a piecewise constant function. 
%   estimaterecovery  - recovery type error estimator.
%   estimaterecovery3 - recovery type error estimator in 3-D.
%   estimateresidual  - residual type error estimator.
%   estimateresidual3 - residual type error estimator in 3-D.
%   faceinterpolate   - interpolate to face elements (RT0 or BDM1).
%   getH1error        - H1 norm of the approximation error.
%   getH1error3       - H1 norm of the approximation error in 3-D.
%   getHcurlerror3NE  - Hcurl norm of approximation error for the lowest order Nedelect element in 3-D.
%   getHcurlerror3NE1 - Hcurl norm of approximation error for the linear Nedelect element in 3-D.
%   getHcurlerror3NE2 - Hcurl norm of approximation error for the quadratic (1st type) Nedelect element in 3-D.
%   getHdiverror3RT0  - Hdiv norm of the approximation error for RT0 in 3-D.
%   getHdiverrorBDM1  - Hdiv norm of the approximation error for BDM1.
%   getHdiverrorRT0   - Hdiv norm of the approximation error for RT0.
%   getL2error        - L2 norm of the approximation error.
%   getL2error3       - L2 norm of the approximation error in 3-D.
%   getL2error3NE     - L2 norm of the lowest order Nedelec edge element.
%   getL2error3NE1    - L2 norm of the linear Nedelec edge element.
%   getL2error3NE2    - L2 norm of the quadratic (1st type) Nedelec edge element.
%   getL2error3RT0    - L2 norm of RT0 element in 3D.
%   getL2errorBDM1    - L2 norm of BDM1 element
%   getL2errorRT0     - L2 norm of RT0 element.
%   gradbasis         - gradient of barycentric basis. 
%   gradbasis3        - gradient of barycentric basis in 3-D.
%   gradu             - gradient of a finite element function.
%   gradu3            - gradient of a finite element function in 3-D.
%   mark              - mark element.
%   nodeinterpolate   - interpolate a piecewise linear function.
%   quadpts           - quadrature points in 2-D.
%   quadpts1          - quadrature points in 1-D.
%   quadpts3          - quadrature points in 3-D.
%   recovery          - recovery a piecewise constant function to a piecewise linear one.
%   recovery3         - recovery a piecewise constant function to a piecewise linear one in 3-D.
%   setboundary       - set type of boundary edges.
%   setboundary3      - set type of boundary faces in 3-D.
%   uniformbisect     - uniformly bisect a 2-D triangulation. 
%   uniformbisect3    - uniformly bisect a 3-D triangulation. 
%   uniformrefine     - uniformly refine a 2-D triangulation.
%   uniformrefine3    - uniformly refine a 3-D triangulation.
%   verifyquadpts     - examples and verfication on quadrature rules.
%   verifyquadpts1    - examples and verfication on quadrature rules in 1-D.
%   verifyquadpts3    - examples and verfication on quadrature rules in 3-D.
