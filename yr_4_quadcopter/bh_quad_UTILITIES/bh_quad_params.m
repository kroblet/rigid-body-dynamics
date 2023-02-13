
%% Propeller props - the properties are for a SINGLE propeller
P_one_prop.LIFT_per_w2   = 0.9e-5;      % (N  /( (rad/sec)^2 ))
P_one_prop.TWIST_per_w2  = 1.2e-7;      % (N.m/( (rad/sec)^2 ))
P_one_prop.thrust_max    = 4.5;        % (N), max thrust that 1 propeller can generate
%P_props.w_max         = 108*2*pi;    % (rad/sec) the max possible propeller rotational speed
P_one_prop.w_max         = sqrt(P_one_prop.thrust_max/P_one_prop.LIFT_per_w2);    % (rad/sec) the max possible propeller rotational speed
P_one_prop.TQ_max        = P_one_prop.TWIST_per_w2*(P_one_prop.w_max.^2); % (N.m) 

%% MOTOR alone params
% Some typical motors used on quadcopters are here:
% http://www.hobbyking.com/hobbyking/store/__21644__18_11_2000kv_Micro_Brushless_Outrunner_10g_AU_Warehouse_.html?strSearch=2000kv
% http://www.hobbyking.com/hobbyking/store/__46010__2728_Brushless_Outrunner_Motor_1000kv_AU_Warehouse_.html?strSearch=1000kv

P_mot_alone.ID = 0;          % (-)           ID
P_mot_alone.Ka = 0.010808;   % (N.m/A)       Torque constant
P_mot_alone.b  = 6.0961E-08; % (N.m.sec/rad) Viscous friction co-efficient
P_mot_alone.I  = 1.6994E-07; % (kg.m^2)      Rotor ONLY Inertia
P_mot_alone.R  = 1.07;       % (Ohms)        Armature resistamce
P_mot_alone.L  = 6.0793E-05; % (Henry)       Armature Inductance
P_mot_alone.Ke = 0.0051;     % (V.sec/rad)   Back EMF constant


%% MOTOR with propeller params
P_mot   = P_mot_alone;
P_mot.I = 0.13198e-4;     % (kg.m^2)      Rotor+Prop Inertia

for kk=1:4
   P_mot_vec(kk,1)    = P_mot;
   P_mot_vec(kk,1).ID = kk;
end

%% Electric supply params
P_elec.V_mot_min = 0;     % (volts)  Min voltage that a motor can receive
P_elec.V_mot_max = 9;     % (volts)  Max voltage that a motor can receive

%% Vehicle params - part 1
% The inertia of the airframe AND propellers
% REF:  bh_do_SYSTEM_m_and_I_calcs_3bladed_prop.mlx 
%       bh_do_ROTOR_I_calcs_3blade.mlx 

P_veh.I = [ ...
   0.005831943165131                   0                   0;
                   0   0.005831943165131                   0;
                   0                   0   0.011188595733333;
           ]; % (kg.m^2)
       
P_veh.mass        = 0.9272;    % (kg)
P_veh.L_CG_to_mot = 0.2;     % (m)
P_veh.Xe_init     = [0;0;0]; % (m)       Initial position in INERTIAl axes
P_veh.Vb_init     = [0;0;0]; % (m/sec)   Initial velocity in BODY axes
P_veh.wb_init     = [0;0;0]; % (rad/sec) Initial body rates
P_veh.eul_init    = [0;0;0]; % (rad)     Initial EULER angles [yaw,pitch,roll]

%% Combined Rotor and prop params
% The inertia of the motor shaft and 3-blade propeller assembly
% REF:  bh_do_ROTOR_I_calcs_3blade.mlx  

P_rot.I =  1.0e-04 * ...
[   0.141353444756554                   0                   0;
                    0   0.141353444756554                   0;
                    0                   0   0.131984333333333;
];       
P_rot.mass        = 17.8/1000; % (kg)
P_rot.Izz         = P_rot.I(3,3);

%% Vehicle airframe params
% The inertia of the airframe.  The propellers are NOT included
% REF:  bh_do_SYSTEM_m_and_I_calcs_3bladed_prop.mlx 
%       bh_do_ROTOR_I_calcs_3blade.mlx 
P_veh_no_props   = P_veh;

P_veh_no_props.I = [ ...
   0.005594125770533       0                   0;
   0                       0.005594125770533   0;
   0                       0                   0.010773333333333; ];

P_veh_no_props.mass        = P_veh.mass - 4*P_rot.mass;
P_veh_no_props.H_CG_to_mot = (0.5+6-1.96)/100; %m

%% Vehicle tilt angle limits:
% compute the maximum pitch and roll angle that we will allow the vehicle
% to make.
P_veh.LIMITS_tilt_max_deg = bh_calc_tilt_max( P_veh.mass, 4*P_one_prop.thrust_max );

%% define some WAYPOINT data
%       O   A    B   C    D   E    F    G   H     I    J    K    L    O
X_WP = [0,  5,  25,  25   5,  5,  25,  25,  5,    5,  25,  25,   0,   0];
Y_WP = [0,  5,   5,  10, 10, 15,  15   20, 20,   25,  25,  30,  30,   0];
T_WP = [0,  7,  27,  32, 52, 57,  77,  82, 102, 107, 127, 132, 157, 187];
Z_WP = [0,  5,   5,   5,  5,  5,   5,   5,   5,   5,   5,   5,   5,   0];

X_TS_OBJ = timeseries(X_WP, T_WP);
Y_TS_OBJ = timeseries(Y_WP, T_WP);
Z_TS_OBJ = timeseries(Z_WP, T_WP);

clear X_WP Y_WP Z_WP T_WP

%% define another set of WAYPOINT data 
% a circle - used for Simscape

XYZ_circ_mat = LOC_define_simple_circle_waypoints();

%% Local helper functions
function XYZ_mat = LOC_define_simple_circle_waypoints()

    R         = 2;
    theta_deg = 0:5:360;

    X = R*cosd(theta_deg);
    Y = R*sind(theta_deg);
    Z = ones(size(X));

    XYZ_mat = [ X(:), Y(:), Z(:) ];

    XYZ_mat = [    0,   0, 0;
                 0.2, 0, 0;
                 0.4, 0, 0;
                 0.6, 0, 0;
                 0.8, 0, 0;                 
                 1.0, 0, 0;
                 1.2, 0, 0;
                 1.4, 0, 0;                 
                 1.6, 0, 0;                 
                 1.8, 0, 0;                 
                 2.0, 0, 0;                 
                 2.0, 0, 0.2;
                 2.0, 0, 0.4;
                 2.0, 0, 0.6;
                 2.0, 0, 0.8;                 
                 %2.0, 0, 1.0;       %_###_BH_07_Nov_2019 - R2019b issue fix                          
                 XYZ_mat ];
end