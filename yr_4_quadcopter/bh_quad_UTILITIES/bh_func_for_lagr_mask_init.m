function [  Iair_xx_yy_zz,  ...
            Iair_yz_zx_xy,  ...
            Iprop_xx_yy_zz, ...
            Iprop_yz_zx_xy, ...
            xme_dot_0,      ...
            eul_dot_0  ] =  ...
                            bh_func_for_lagr_mask_init( air_inertia_0,  ... 
                                                        prop_Imat,      ...
                                                        Vm_0,           ...
                                                        pm_0,           ...
                                                        eul_0)
% NOTE:
%  Vm_0  : initial velocity in body fixed components
%  pm_0  : initial angular rates in body fixed components
%  eul_0 : initial Euler angles

Iair_xx_yy_zz = [air_inertia_0(1,1), ...
                 air_inertia_0(2,2), ...
                 air_inertia_0(3,3) ];
             
Iair_yz_zx_xy = [air_inertia_0(2,3), ...
                 air_inertia_0(3,1), ...
                 air_inertia_0(1,2) ];  
             
             
Iprop_xx_yy_zz = [prop_Imat(1,1), ...             
                  prop_Imat(2,2), ...
                  prop_Imat(3,3)  ];
              
Iprop_yz_zx_xy = [prop_Imat(2,3), ...             
                  prop_Imat(3,1), ...
                  prop_Imat(1,2)  ]; 
%-------------------------------------
phi   = eul_0(1);
theta = eul_0(2);
psi   = eul_0(3);

p     = pm_0(1);
q     = pm_0(2);
r     = pm_0(3);

bRg = [ ...
                              cos(phi)*cos(theta),                              cos(theta)*sin(phi),         -sin(theta);
 cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(theta)*sin(psi);
 sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), cos(psi)*cos(theta);
 ];

euler_rates = [(r.*cos(psi)+q.*sin(psi))./cos(theta);
                q.*cos(psi)-r.*sin(psi);
               (p.*cos(theta)+r.*cos(psi).*sin(theta)+q.*sin(psi).*sin(theta))./cos(theta)];
           
           
xme_dot_0 = (bRg .') * Vm_0(:);
eul_dot_0 = euler_rates;
end

