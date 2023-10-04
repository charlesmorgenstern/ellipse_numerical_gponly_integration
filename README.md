# ellipse_numerical_gponly_integration
Numerical methods for GP only integration of the ellipse system


Functions:

getpaths(n,a,h) generate gradient paths. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

plotpaths(n,a,h) plot level sets and gradient paths. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

curv(r) get curvature of level set at point r with exact formula

curvfdm2(r,h) get curvature of level set at point r with finite differences. h is step size.

getpathlengths(gps) calculate the length of gradient path segments with Euclidean distance. gps are gradient paths.

getpathlengthsex(gps) calculate the length of gradient path segments with integral. gps are gradient paths.

checkpathlengths(n,a,h) calculate maximum relative error of gradient path segment lengths. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getcurvs(gps) get curvature of level sets at the needed points using finite difference curvature code. gps are gradient paths.

getcurvsex(gps) get curvature of level sets at the needed points using exact curvature code. gps are gradient paths.

getcurvderivatives(gps) get derivative of level set curvature in direction of level set with finite difference. gps are gradient paths.

getcurvderivativesex(gps) get derivative of level set curvature in direction of level set with exact formula. gps are gradient paths.

getpathcurv(gps) get curvature of gradient paths with finite differences. gps are gradient paths.

getpathcurvex(gps) get curvature of gradient paths with exact formula. gps are gradient paths.

getleftpathcurv(gps) get gradient path curvature of neighboring gradient path at corresponding points. gps are gradient paths.

getleftpathlengths(gps) get lengths of corresponding gradient path segments of neighboring gradient path. gps are gradient paths.

getarclengths(gps,curvs,lengths,lengths2,dk,a,k,k2) get arc lengths of level set segments. gps are gradient paths. curvs are level set curvatures. lengths are gradient path segment lengths. lengths2 are gradient path segment lengths of neighboring gradient path. a is x-intercept of staring ellipse. dk are the derivatives of the level set curvatures. k are gradient path curvatures. k2 are gradient path curvatures of neighboring gradient path.

getarea(n,a,h) get area of ellipse using gradient path only integration and all numerical methods.  n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getareagb(n,a,h) get area of individual gradient bundles with gradient path only integration and all numerical methods. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getareagreens(n,a,h) get area of individual gradient bundles using Green's theorem. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getareaerror(n,a,h) get relative error of ellipse area using gradient path only integration and all numerical methods. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getn(gps) get maximum number of points along a gradient path. gps are gradient paths.

errortablearea(a,nn,nstart,hstart) generate table of relative error and convergence of ellipse area as number of gradient paths is increased. a is x-intercept of starting ellipse. nn is number of refinements. nstart is number of gradient paths for coarsest case. hstart is step size for gradient path generation and integration for coarsest case. 

plotarclengtherror(n,a,h) plot relative error of level set arc lengths on largest ellipse. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

animatearclengtherror() animate plots of arc length error as starting ellipse size changes.

animatesetup() animate gradient paths as starting ellipse size changes.

getc(gps) get parameters c of parabolas y=cx^2 which are gradient paths. gps are gradient paths.

getgperrorL2(gps) get L2 error of gradient paths to check RK4. gps are gradient paths.

plotpathswexact(n,a,h) plot RK4 gradient paths and exact gradient paths on same plot

getarclengthsex(gps) get exact arc lengths of level set segments. gps are gradient paths.

getarclengtherrors(n,a,h) get relative error of arc lengths of level set segments. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

plotarclengtherrors(n,a,h) plot relative error of the arc lengths of level set segments. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

plotgberror(n,a,h) plot relative error of individual gradient bundle areas title with number of gradient bundles. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

animategberror() animate plots of relative error of individual gradient bundle areas with varying number of gradient bundles.

plotgberror2(n,a,h) plot relative error of individual gradient bundle areas title with starting ellipse size. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

animategberror2() animate plots of relative error of individual gradient bundle areas with varying starting ellipse size.

plotcurverrors(n,a,h) plot relative errors of level set curvatures. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getcurverrors(n,a,h) get relative error of level set curvatures. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getcurvsex(gps) get exact level set curvatures. gps are gradient paths.

getpatherrors(n,a,h) get relative error of gradient path segment lengths. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

plotpatherrors(n,a,h) plot relative errors of gradient path segment lengths. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getpathcurverrors(n,a,h) get realtive errors of gradient path curvatures. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

plotpathcurverrors(n,a,h) plot relative errors of gradient path curvatures. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

getcurvderivativeerrors(n,a,h) get relative errors of level set curvature derivatives. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.

plotcurvderivativeerrors(n,a,h) plot relative errors of level set curvature derivatives. n is number of gradient paths. a is x-intercept of starting ellipse. h is step size for path generation.


