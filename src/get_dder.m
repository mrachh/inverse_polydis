function [dder] = get_dder(zk,vs,nvs,vgrad,angs,targs,h,u)

    vnrm = norm(vgrad,'fro');
    vts = vs;
    vts = vts + h*vgrad/vnrm;
       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
    vp = norm(u-u_s,'fro')^2;
    
    vts = vs;
    vts = vts - h*vgrad/vnrm;
       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
    vm = norm(u-u_s,'fro')^2;
    
    vts = vs;
       xy_s = sum(vts')'/nvs;
 [u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vts,zk,targs,angs,xy_s);
    vmid = norm(u-u_s,'fro')^2;
    dder = (vm + vp - 2*vmid)/(h*h);
end

