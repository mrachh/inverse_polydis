function [u,varargout] = helm_dirichlet_solver(verts,zk,targs,angs,xyin,cparams)
%HELM_DIRICHLET_SOLVER solves the helmholtz Dirichlet value problem
%for an analytic solution and a collection of plane waves for a rounded
%polygon
%
% Syntax: [u,varargout] = helm_dirichlet_solver(verts,zk,targs,angs,xyin,cparams)
% Input:
%   verts = xy coordinates of vertices describing the polygon (2,nverts)
%   zk - Helmholtz wave number
%   targs - xy coordinates of targets in exterior for evaluating potential
%   (2,nt)
%   angs - angles of incidence for the incoming plane waves (nangs,1)
%   xyin - location of interior point 
% Optional input:
%   cparams - options structure
%       cparams.rounded = true if corner rounding is to be used.
%                         false if no rounding is used (true)
%       cparams.autowidths = automatically compute widths (true)
%       cparams.autowidthsfac = if using autowidths, set widths
%                             to autowidthsfac*minimum of adjoining
%                             edges (0.4)
% 	    cparams.eps - resolve curve to tolerance eps
%                    resolve coordinates, arclength,
%          	     and first and second derivs of coordinates
%		         to this tolerance (1.0e-6)
% output:
%   u - potential at targets (nt,nangs)
%   optional output arguments:
%   chnkr - discretized chunker used in solver
%   bd_sol - solution of integral equation corresponding to incident plane
%       waves (n,nangs)
%   F - compressed linear system corresponding to integral equation using
%       rskelf
%   err - estimated error in solving dirichlet problem
%

if(nargin<=5)
    cparams_use = [];
    cparams_use.rounded = true;
    cparams_use.autowidths = true;
    cparams_use.autowidthsfac = 0.4;
    cparams_use.eps = 1e-6;
else
    cparams_use = cparams;
end
pref.k = 16;
pref.dim = 2;
chnkr = chunkerpoly(verts,cparams_use,pref);
    
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'C',1);
dval = 0.5;
opts_flam = [];
opts_flam.flamtype = 'rskelf';
opts_flam.rank_or_tol = 1e-10;
F = chunkerflam(chnkr,fkern,dval,opts_flam);

% Test exact solution
srcinfo = []; srcinfo.r = xyin; targinfo = []; targinfo.r = chnkr.r;
targinfo.r = reshape(targinfo.r,2,chnkr.k*chnkr.nch);

bd_test = chnk.helm2d.kern(zk,srcinfo,targinfo,'s');

bd_sol_ex = rskelf_sv(F,bd_test);
varargout{1} = chnkr;
varargout{3} = F;
utest = chunkerkerneval(chnkr,fkern,bd_sol_ex,targs);

targinfo = [];
targinfo.r = targs;
uex = chnk.helm2d.kern(zk,srcinfo,targinfo,'s');

err = norm(utest-uex,'fro')/norm(uex,'fro');
varargout{4} = err;


rval = chnkr.r;
rval = reshape(rval,2,chnkr.k*chnkr.nch);
xval = rval(1,:);
yval = rval(2,:);
uinc = -exp(1i*zk*(xval'*cos(angs)' + yval'*sin(angs)'));
bd_sol = rskelf_sv(F,uinc);
varargout{2} = bd_sol;
nt = length(targs);
nangs = length(angs);
u = complex(zeros(nt,nangs));
for i=1:nangs
    u(:,i) = chunkerkerneval(chnkr,fkern,bd_sol(:,i),targs);
end

   
end