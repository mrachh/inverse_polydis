function [vs,nvs_out] = add_verts(vs_in,nvs,zk,u,targs,angs,nangs)

inds = [];
vs = vs_in;
nadd = 0;
 for i=1:nvs-1
     rad = norm(vs(:,i+1+nadd)-vs(:,i+nadd),'fro');
     if (rad >2*2*pi/zk)
        disp('here')
        pmid =  (vs(:,i+1+nadd)+vs(:,i+nadd))/2;
        vs = [vs(:,1:(i+nadd)),pmid,vs(:,(i+1+nadd):end)];
        inds = [inds,i+nadd+1];
        nadd = nadd + 1;
     end
 end

rad = norm(vs(:,end)-vs(:,1),'fro');
if (rad >2*2*pi/zk)
   disp('here')
   pmid =  (vs(:,end)+vs(:,1))/2;
   vs = [vs(:,1:end),pmid];
   nadd = nadd + 1;
   inds = [inds,nvs+nadd];
end

nvs_out = nvs + nadd;

vs
nvs_out

xy_s = sum(vs')'/nvs_out;
[u_s,u_grad,chnkr_s,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(vs,zk,targs,angs,xy_s);

ugrad_mat = reshape(u_grad,[nangs*nangs,2*nvs_out]);
vdelt     = reshape(u_s-u,[nangs*nangs,1]);
vgrad = reshape(2*real(vdelt'*ugrad_mat),[2,nvs_out]);
ugrads = sum(abs(vgrad).^2,1);
ugrads = sqrt(ugrads)

if (numel(inds) >0)
inds
[m,i] = max(ugrads(inds))
inds(i(1)) = []
vs(:,inds) = []
nvs_out = nvs + 1;
else
vs = vs_in;
nvs_out = nvs;
end    


