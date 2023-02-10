classdef bh_4bar_kin_circ_CLS
%_#########################################################################
% EXPECTED USAGE:
%    obj               = bh_4bar_kin_circ_CLS(L1,L2,L3,L4,theta_rad_SCALAR)
% 
%    tf_is_valid_pose  = obj.is_valid_pose();
%    
%    [ phi_pair_rad, ...
%      alpha_pair_rad] = obj.calc_phi_alpha_rad();   
%_#########################################################################
    properties (SetAccess = protected)
        L1        = 0;
        L2        = 0;
        L3        = 0;
        L4        = 0;
        theta_rad = 0;
    end
    
    methods
        function obj = bh_4bar_kin_circ_CLS(L1,L2,L3,L4,theta_rad)
          % ALLOWED USAGE:
          %   >> obj = bh_4bar_kin_circ_CLS(L1,L2,L3,L4,theta_rad)

             if(5~=nargin)
                error("###_ERROR:  we need EXACTLY 5 inputs --->  L1, L2, L3, L4, theta_rad");
             end
          
             % *ONLY* consider SCALR inputs that are FINITE and NOT NaN
             validateattributes(L1,        {'double'},{'scalar', 'finite','nonnan'});
             validateattributes(L2,        {'double'},{'scalar', 'finite','nonnan'});
             validateattributes(L3,        {'double'},{'scalar', 'finite','nonnan'});
             validateattributes(L4,        {'double'},{'scalar', 'finite','nonnan'});
             validateattributes(theta_rad, {'double'},{'scalar', 'finite','nonnan'});
             
             % OK - we're good to go
             obj.L1        = L1;
             obj.L2        = L2;
             obj.L3        = L3;
             obj.L4        = L4;  
             obj.theta_rad = theta_rad;
        end
        %------------------------------------------------------------------
        function tf_is_valid_pose = is_valid_pose(obj)
            [xb_col, yb_col ]  = obj.calc_ORIG_xb_col_AND_yb_col();            
            [~, ~, tf_struct ] = obj.clean_xB_and_yB_vectors(xb_col,yb_col);           
            tf_is_valid_pose   = tf_struct.tf_is_valid; 
            
            assert( numel(tf_is_valid_pose)==1 )
        end
        %------------------------------------------------------------------
        function [xb_CLEAN_col, yb_CLEAN_col] = calc_CLEAN_xb_col_AND_yb_col(obj)
        
            % *** ATTENTION ***
            %   the designed behaviour is that an INVALID machine will 
            %   return NaN values for the angles.
            
            [xb_col, yb_col]  = obj.calc_ORIG_xb_col_AND_yb_col();
            
            [xb_CLEAN_col, ...
             yb_CLEAN_col, ...
             tf_struct      ] = obj.clean_xB_and_yB_vectors(xb_col,yb_col);
        end
        %------------------------------------------------------------------
        function [phi_pair_rad, alpha_pair_rad] = calc_phi_alpha_rad(obj) 
            
            [xb_CLEAN_col, yb_CLEAN_col] = obj.calc_CLEAN_xb_col_AND_yb_col();
            
            assert( numel(xb_CLEAN_col)==2 )
            assert( numel(yb_CLEAN_col)==2 )
            
            [L1,L2,L3,L4]   = obj.get_L();
            theta_deg       = obj.get_theta_DEG();            
            L_val_vec       = [L1,L2,L3,L4];
            
            [phi_1_deg_val, alpha_1_deg_val] = ...
                   obj.calc_phi_alpha_deg_UNIT(xb_CLEAN_col(1), ...
                                               yb_CLEAN_col(1), ...
                                               theta_deg, L_val_vec);
                                           
            [phi_2_deg_val, alpha_2_deg_val] = ...
                   obj.calc_phi_alpha_deg_UNIT(xb_CLEAN_col(2), ...
                                               yb_CLEAN_col(2), ...
                                               theta_deg, L_val_vec);

             phi_pair_rad   = deg2rad([  phi_1_deg_val;   phi_2_deg_val]);                              
             alpha_pair_rad = deg2rad([alpha_1_deg_val; alpha_2_deg_val]);
        end
        %------------------------------------------------------------------
        function [xb_col, yb_col] = calc_ORIG_xb_col_AND_yb_col(obj)
            
            [L1, L2, L3, L4] = obj.get_L();
            theta_deg        = obj.get_theta_DEG();
            [xA,yA]          = obj.get_xA_yA();

            xb_col = LOC_mlf_xB(L1,L2,L3,L4,xA,yA);
            yb_col = LOC_mlf_yB(L1,L2,L3,L4,xA,yA);
            
            xb_col = xb_col(:);
            yb_col = yb_col(:);
            
            assert( numel(xb_col)==2 )
            assert( numel(yb_col)==2 )
        end
        %------------------------------------------------------------------        
    end % METHODS
%_#########################################################################
   methods(Access = protected)
       function [L1,L2,L3,L4] = get_L(obj)
            L1 = obj.L1;
            L2 = obj.L2;
            L3 = obj.L3;
            L4 = obj.L4;
       end
       %-------------------------------------------------------------------
       function theta_deg = get_theta_DEG(obj)
                theta_deg = rad2deg(obj.theta_rad);
       end
       %-------------------------------------------------------------------
       function [xA,yA] = get_xA_yA(obj)
           [L1,~,~,~]   = get_L(obj);
            theta_deg   = get_theta_DEG(obj);            
            xA          = L1*cosd(theta_deg);
            yA          = L1*sind(theta_deg);
       end
       %-------------------------------------------------------------------
        function [xB_vec_real, yB_vec_real, tf_struct] = clean_xB_and_yB_vectors(obj, the_xB_vec,  the_yB_vec)
            
            % *** ATTENTION ***
            %   the designed behaviour is that an INVALID machine will 
            %   return NaN values for the angles.

            validateattributes(the_xB_vec,{'double'},{'vector', 'numel',2});  
            validateattributes(the_yB_vec,{'double'},{'vector', 'numel',2});   

            tf_xB_is_complex_vec  = abs(  imag(the_xB_vec)  ) > 1e-10;
            tf_yB_is_complex_vec  = abs(  imag(the_yB_vec)  ) > 1e-10;  

            tf_is_complex_SCAL    = tf_xB_is_complex_vec(1) | tf_xB_is_complex_vec(2) | ...
                                    tf_yB_is_complex_vec(1) | tf_yB_is_complex_vec(2) ;

            tf_is_real_SCAL       = ~tf_is_complex_SCAL;   

            % So let's create some new matrices with NaN at  
            % the locations where we had INVALID values for xB and yB

            if(tf_is_real_SCAL==true)
                    xB_vec_real = the_xB_vec;
                    yB_vec_real = the_yB_vec;    
            else
                    xB_vec_real = NaN(size(the_xB_vec));
                    yB_vec_real = NaN(size(the_yB_vec));
            end

            % And let's also test for the presence of Inf
            tf_xB_is_inf_vec = isinf(the_xB_vec);
            tf_yB_is_inf_vec = isinf(the_yB_vec);
            tf_is_inf_SCAL   = tf_xB_is_inf_vec(1) | tf_xB_is_inf_vec(2) | ...
                               tf_yB_is_inf_vec(1) | tf_yB_is_inf_vec(2) ;

            % So clean our matrices:
            if(tf_is_inf_SCAL==false)
                    xB_vec_real = the_xB_vec;
                    yB_vec_real = the_yB_vec;    
            else
                    xB_vec_real = NaN(size(the_xB_vec));
                    yB_vec_real = NaN(size(the_yB_vec));
            end    

            % So let's define the notion of a valid pose 
            tf_is_valid_SCAL = tf_is_real_SCAL & ~tf_is_inf_SCAL;       

            % take care of the tf struct output
            tf_struct.tf_is_real  = tf_is_real_SCAL;
            tf_struct.tf_is_inf   = tf_is_inf_SCAL;
            tf_struct.tf_is_valid = tf_is_valid_SCAL;
        end       
       %-------------------------------------------------------------------
       function [phi_deg_val, alpha_deg_val] = ...
                   calc_phi_alpha_deg_UNIT(obj, xB, yB , theta_deg, L_val_vec)

              validateattributes(xB,       {'double'},{'scalar'});
              validateattributes(yB,       {'double'},{'scalar'});
              validateattributes(theta_deg,{'double'},{'scalar'});  
              validateattributes(L_val_vec,{'double'},{'vector', 'numel',4});  

              phi_deg_val   = NaN;
              alpha_deg_val = NaN;

              if( ~isfinite(xB) | ~isfinite(yB) | ~isreal(xB) | ~isreal(yB) | ~isfinite(theta_deg) )
                 return
              end

              % so if we're here ... we can proceed with a calculation
              L1 = L_val_vec(1);
              L2 = L_val_vec(2);
              L3 = L_val_vec(3);
              L4 = L_val_vec(4);  

              xA = L1*cosd(theta_deg);
              yA = L1*sind(theta_deg);

              xD = L4;
              yD = 0;

              % calc phi_deg 
              vec_D2B       = [xB;yB] - [xD;yD];
              phi_deg_val = atan2d(vec_D2B(2), vec_D2B(1) );

              % calc alpha_deg  
              vec_A2B       = [xB;yB] - [xA;yA];  
              alpha_deg_val = atan2d(vec_A2B(2), vec_A2B(1) );  
        end       
       %-------------------------------------------------------------------
   end % methods(Access = protected)
end % CLASSDEF
%_#########################################################################
% Support functions
%_#########################################################################
function xB_sol = LOC_mlf_xB(L1,L2,L3,L4,xA,yA)
%TMP_MLF_XB
%    XB_SYM_SOL = TMP_MLF_XB(L1,L2,L3,L4,XA,YA)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    26-Sep-2019 09:13:59

t2 = L2.^2;
t3 = L3.^2;
t4 = L4.^2;
t5 = L4.^3;
t6 = xA.^2;
t7 = xA.^3;
t8 = yA.^2;
t9 = L4.*xA.*2.0;
t10 = L2.*L3.*2.0;
t11 = -t9;
t12 = L4.*t2;
t13 = L4.*t3;
t14 = t2.*xA;
t15 = t3.*xA;
t16 = L4.*t6;
t17 = t4.*xA;
t18 = L4.*t8;
t19 = t8.*xA;
t20 = -t2;
t21 = -t3;
t22 = -t4;
t23 = -t6;
t24 = -t8;
t25 = -t13;
t26 = -t14;
t27 = -t16;
t28 = -t17;
t29 = t4+t6+t8+t11;
t32 = t2+t3+t9+t10+t22+t23+t24;
t30 = 1.0./t29;
t31 = t10+t20+t21+t29;
t33 = t31.*t32;
t34 = sqrt(t33);
t35 = t34.*yA;
xB_sol = [(t30.*(t5+t7+t12+t15+t18+t19+t25+t26+t27+t28-t35))./2.0;(t30.*(t5+t7+t12+t15+t18+t19+t25+t26+t27+t28+t35))./2.0];

end
%--------------------------------------------------------------------------
function yB_sol = LOC_mlf_yB(L1,L2,L3,L4,xA,yA)
%TMP_MLF_YB
%    YB_SYM_SOL = TMP_MLF_YB(L1,L2,L3,L4,XA,YA)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    26-Sep-2019 09:14:00

t2 = L2.^2;
t3 = L3.^2;
t4 = L4.^2;
t5 = xA.^2;
t6 = yA.^2;
t7 = yA.^3;
t8 = L4.*xA.*2.0;
t9 = L2.*L3.*2.0;
t21 = L4.*xA.*yA.*-2.0;
t10 = -t8;
t11 = t2.*yA;
t12 = t3.*yA;
t13 = t4.*yA;
t14 = t8.*yA;
t15 = t5.*yA;
t16 = -t2;
t17 = -t3;
t18 = -t4;
t19 = -t5;
t20 = -t6;
t22 = -t11;
t23 = t4+t5+t6+t10;
t26 = t2+t3+t8+t9+t18+t19+t20;
t24 = 1.0./t23;
t25 = t9+t16+t17+t23;
t27 = t25.*t26;
t28 = sqrt(t27);
t29 = L4.*t28;
t30 = t28.*xA;
yB_sol = [(t24.*(t7+t12+t13+t15+t21+t22-t29+t30))./2.0;(t24.*(t7+t12+t13+t15+t21+t22+t29-t30))./2.0];
end
%--------------------------------------------------------------------------
