function [hess] = get_hess(zk,vs,nvs,angs,targs,u,h)
    
    hess = zeros([size(vs),size(vs)]);
    for i1=1:nvs
    for j1=1:2
    for i2=i1:nvs
    for j2=1:2    
       vts = vs;
       vts(j1,i1) = vts(j1,i1) + h;
       vts(j2,i2) = vts(j2,i2) + h;
       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
        vpp = norm(u-u_s,'fro')^2;
        
       vts = vs;
       vts(j1,i1) = vts(j1,i1) + h;
       vts(j2,i2) = vts(j2,i2) - h;
       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
        vpm = norm(u-u_s,'fro')^2; 
        
       vts = vs;
       vts(j1,i1) = vts(j1,i1) - h;
       vts(j2,i2) = vts(j2,i2) + h;
       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
        vmp = norm(u-u_s,'fro')^2; 
        
       vts = vs;
       vts(j1,i1) = vts(j1,i1) - h;
       vts(j2,i2) = vts(j2,i2) - h;
       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
        vmm = norm(u-u_s,'fro')^2; 
        
        hess(j1,i1,j2,i2) = (vpp-vmp-vpm+vmm)/(4*h*h);
        hess(j2,i2,j1,i1) = (vpp-vmp-vpm+vmm)/(4*h*h);
    end
    end
    end
    end
end

