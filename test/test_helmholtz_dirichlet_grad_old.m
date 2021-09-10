addpath('../src');
verts = [0,1,1,0;0,0,1,1];
xyin = [0.3;0.4];
zk = 1.1;
angs = [pi/3];
targs = [3.1, 4.2, 4; -1.7,2.5,-3.6];
[u,chnkr,bd_sol,F,err1] = helm_dirichlet_solver(verts,zk,targs,angs,xyin);


% First test evaluation of normal derivative on the boundary
srcinfo = []; srcinfo.r = xyin; targinfo = []; targinfo.r = chnkr.r;
targinfo.r = reshape(targinfo.r,2,chnkr.k*chnkr.nch);
targinfo.d = chnkr.d;
targinfo.d = reshape(targinfo.d,2,chnkr.k*chnkr.nch);
ubdry = chnk.helm2d.kern(zk,srcinfo,targinfo,'s');
dudnbdry = chnk.helm2d.kern(zk,srcinfo,targinfo,'sprime');
sigma = rskelf_sv(F,ubdry);
    
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'s',1);
dval = 0.0;
opts_flam = [];
opts_flam.flamtype = 'rskelf';
opts_flam.rank_or_tol = 1e-10;
F_s = chunkerflam(chnkr,fkern,dval,opts_flam);


fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'sprime',1);
dval = -0.5;
F_sprime = chunkerflam(chnkr,fkern,dval,opts_flam);
 
dp = helm_dprime(zk,chnkr,F_s,sigma);
sp = rskelf_mv(F_sprime,sigma);
 
dudn = dp + 1i*sp;
err_dudn = norm(dudn-dudnbdry);
fprintf('Error in computing dudn corresponding to represenation: %5.2e\n',...
    err_dudn);

% Now start gradient computation
grad = vertgrad(chnkr);
[~,nv] = size(verts);
pert = randn(2,nv);
pert = zeros(2,nv);
pert(1) = 1;
g = grad*pert(:);
gx = g(1:2:end); gy = g(2:2:end);

rn = normals(chnkr);
rn = reshape(rn,[2,chnkr.k*chnkr.nch]);
rnx = rn(1,:)'; rny= rn(2,:)';

% compute dudn corresponding to plane wave data
% currently assumes there is only one plane wave
dp = helm_dprime(zk,chnkr,F_s,bd_sol);
sp = rskelf_mv(F_sprime,bd_sol);

dudn = dp + 1i*sp;

% add in the incident field
rval = chnkr.r;
rval = reshape(rval,2,chnkr.k*chnkr.nch);
xval = rval(1,:);
yval = rval(2,:);
uinc = exp(1i*zk*(xval'*cos(angs)' + yval'*sin(angs)'));
dudninc = uinc.*1i.*zk.*(rnx.*cos(angs)' + rny.*sin(angs)');

dsdt = sqrt(chnkr.d(1,:,:).^2 + chnkr.d(2,:,:).^2);
dsdt = reshape(dsdt, [chnkr.k*chnkr.nch,1]);
% Compute boundary data  = -(gx, gy) \cdot (rnx,rny)*dudn;
bd_data_grad = -(gx.*rnx + gy.*rny).*(dudn+dudninc);

% Solve the integral equation with data bd_data_grad
sigma_grad = rskelf_sv(F,bd_data_grad);

% Now evaluate solution at targets    
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'C',1);
[~,nt] = size(targs);
nangs = length(angs);
u_grad = complex(zeros(nt,nangs));
for i=1:nangs
    u_grad(:,i) = chunkerkerneval(chnkr,fkern,sigma_grad(:,i),targs);
end

% Now compute gradient through finite differences


ng = 3;
u_grad_fd = complex(zeros(nt,nangs,ng));
for i = 1:ng
    eps = 10^(-(i));
    
    verts2 = verts + pert*eps;
    [u2,chnkr2,bd_sol2,F2,err2] = helm_dirichlet_solver(verts2,zk,targs,...
       angs,xyin);
    u_grad_fd(:,:,i) = (u2 - u)/eps;
     
end


[u,u_grad,chnkr,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(verts,zk,targs,angs,xyin);

pert = randn(2,nv);
[n,nangs,~,~] = size(u_grad);
u_grad_test = reshape(u_grad,[n,nangs,2*nv]);
u_grad_test = reshape(u_grad_test, [n*nangs,2*nv]);
u_grad_test = u_grad_test*pert(:);
u_grad_test = reshape(u_grad_test,[n,nangs]);

ng = 3;
u_grad_fd2 = complex(zeros(nt,nangs,ng));
for i = 1:ng
    eps = 10^(-(i));
    
    verts2 = verts + pert*eps;
    [u2,chnkr2,bd_sol2,F2,err2] = helm_dirichlet_solver(verts2,zk,targs,...
       angs,xyin);
    u_grad_fd2(:,:,i) = (u2 - u)/eps;
     
end

