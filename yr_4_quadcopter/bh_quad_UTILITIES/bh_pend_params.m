
%% Pendulum parameters
P_pend.mass   = 0.25;                % kg - mass
P_pend.L_xyz  = [0.02, 0.02, 0.5 ]'; % m  - side lengths
P_pend.T_xyz  = [0,    0,    0.05]'; % m  - position of spherical joint 
                                     %      relative to airframe CoM, in
                                     %      components of the body fixed 
                                     %      airframe axes
                                    
% Inertia of pendulum through its CoM                                    
P_pend.I      = LOC_calc_I_pend(P_pend.mass, P_pend.L_xyz);  % kg.m^2

% Local Euler angles of {P}-frame relative to the {V}-frame
%   - rotation sequence of Z-Y-X
%   - phi, theta, psi
P_pend.eul_loc_vec = (pi/180)*[0,10,10]'; % rad

% END OF SCRIPT
%--------------------------------------------------------------------------
% Start of local functions
%--------------------------------------------------------------------------
function I_pend = LOC_calc_I_pend(m, L_xyz)

Lx = L_xyz(1);
Ly = L_xyz(2);
Lz = L_xyz(3);

I_pend = ...
      m * [  (Ly^2 + Lz^2)/12,    0,                    0; ...
              0,                 (Lx^2 + Lz^2)/12,      0; ...
              0,                  0,                   (Lx^2 + Ly^2)/12;] ; 
                
end
%--------------------------------------------------------------------------
