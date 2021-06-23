function [vgrad] = get_grad_faster(zk,vs,nvs,angs,targs,u,h)
        vts = vs;

       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
        vm = norm(u-u_s,'fro')^2; 
        
    vgrad = zeros(size(vs));
    for i=1:nvs
    for j=1:2    
       vts = vs;
       vts(j,i) = vts(j,i) + h;
       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
        vp = norm(u-u_s,'fro')^2;
 
        vgrad(j,i) = (vp-vm)/(h);
    end
    end
end

