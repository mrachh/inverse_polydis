addpath('../src');
clear
verts = [0,1,1,0;0,0,1,1];
xyin = [0.3;0.4];
zk = 1.1;
angs = [pi/3];
targs = [3.1, 4.2, 4; -1.7,2.5,-3.6];

[u,u_grad,chnkr,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(verts,zk,targs,angs,xyin);
[~,nv] = size(verts);
nangs = length(angs);
n = chnkr.k*chnkr.nch;
[~,nt] = size(targs);
pert = randn(size(verts));
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

err = norm(u_grad_test-u_grad_fd2(:,:,ng))/norm(u_grad_fd2(:,:,2));
fprintf('error in gradient evaluation %5.2e\n',err);