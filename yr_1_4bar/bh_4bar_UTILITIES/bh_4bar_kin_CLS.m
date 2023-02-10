classdef bh_4bar_kin_CLS
%_#########################################################################
% EXPECTED USAGE:
%
%  Initialization:
%    >> OBJ_kin = bh_4bar_kin_CLS(L1, L2, L3, L4)
%    >> OBJ_kin = bh_4bar_kin_CLS(L1, L2, L3, L4, theta_rad)
%
%  VALID machine and pose
%    >> OBJ_kin.is_valid_machine()
%    >> OBJ_kin.is_crank_machine()
%    >> OBJ_kin.is_valid_pose(theta_rad_vec)
%
%  Solution pairs for Z:
%    >> [z_pair, num_pair, den_pair] = OBJ_kin.calc_z(theta_rad)
%
%  Solution pairs for PHI and ALPHA
%    >> [phi_pair_rad, alpha_pair_rad] = OBJ_kin.calc_phi_alpha_rad(theta_rad)
%
%  Positive PHI pose:
%    >> OBJ_kin.plot_phi_pos();
%    >> [phi_rad, alpha_rad] = OBJ_kin.get_phi_alpha(theta_rad);
%
%  ROOTS for theta sweep:
%    >> [z_pair_mat, theta_deg_col, tf_is_real] = OBJ_kin.calc_z_values();
%
%  Visualizations:
%    >> OBJ_kin.plot()
%    >> OBJ_kin.plot_angles_flat()
%    >> OBJ_kin.plot_angles_flat()
%    >> OBJ_kin.plot_angles_round()
%    >> OBJ_kin.plot_z()
%_#########################################################################
    properties (SetAccess = protected)
        L1 = 0;
        L2 = 0;
        L3 = 0;
        L4 = 0;
        theta_rad = 0;
    end
    %----------------------------------------------------------------------
    properties (SetAccess = protected)
      %tgt_theta_deg_col    = [(-180:1:0), (1:180)]';
      tgt_theta_deg_col    = [(-180:0.5:0), (0.5:0.5:180)]';
      tf_tgt_is_valid_pose_col = false(size(361,1)); 
    end
    %----------------------------------------------------------------------
%     properties (SetAccess = protected)
%        sweep_stuff;
%     end
    %---------------------------------------------------------------------- 
    methods
        function obj = bh_4bar_kin_CLS(L1,L2,L3,L4,theta_rad)
          % ALLOWED USAGE:
          %   >> obj = bh_4bar_kin_CLS()
          %   >> obj = bh_4bar_kin_CLS(L1,L2,L3,L4)
          %   >> obj = bh_4bar_kin_CLS(L1,L2,L3,L4,theta_rad)

             if(0==nargin)
                L1 = 0.4;
                L2 = 1.4;
                L3 = 1.2;
                L4 = 1;
             end
          
             obj.L1 = L1;
             obj.L2 = L2;
             obj.L3 = L3;
             obj.L4 = L4;

             % compute the valid pose angles for theta
             theta_rad_vec                = deg2rad(obj.tgt_theta_deg_col);
             obj.tf_tgt_is_valid_pose_col = obj.is_valid_pose(theta_rad_vec);

             % make sure we have some valid angles for theta
             if(all(~obj.tf_tgt_is_valid_pose_col))
                 obj.plot_circles();
                 error('###_ERROR:  this machine has NO valid poses !');
             else
                 % find the first valid pose
                 tmp_k                 = find(obj.tf_tgt_is_valid_pose_col, 1);
                 first_valid_theta_rad = deg2rad(obj.tgt_theta_deg_col(tmp_k)); 
             end

             % make sure that the initial theta value is valid
             if(nargin < 5)
                 obj.theta_rad = first_valid_theta_rad;
             else
                 % we supplied theta_rad as an input
                 obj.theta_rad = theta_rad;
             end

             % check the pose of the machine
             if( ~obj.is_valid_pose(obj.theta_rad) )
                 obj.plot_valid_theta();
                 obj.plot_angles_flat();
                 obj.plot_circles();
                error('###_ERROR:  your specified theta is invalid !') 
             end
        end % bh_4bar_kin_CLS()
        %------------------------------------------------------------------
        function tf_val = is_valid_machine(obj)
           % To be a valid machine:
           %  a.) there must be *SOME* values of theta for which the 
           %      machine's pose can be computed. 
            
           if( any(obj.tf_tgt_is_valid_pose_col) )  
               tf_val = true;
           else
               tf_val = false;
           end
        end
        %------------------------------------------------------------------
        function tf_val = is_crank_machine(obj)
           % To be a CRANK machine:
           %  a.) ALL values of theta should produce a valid machine pose 
            
           if( all(obj.tf_tgt_is_valid_pose_col) )  
               tf_val = true;
           else
               tf_val = false;
           end
        end
        %------------------------------------------------------------------
        function tf_is_valid_pose_col = is_valid_pose(obj, theta_rad_vec)
           N               = length(theta_rad_vec); 
           tf_is_valid_col = false(N,1);
           
           L_1 = obj.L1;
           L_2 = obj.L2;
           L_3 = obj.L3;
           L_4 = obj.L4;
            
          for kk=1:length(theta_rad_vec)
              theta_rad      = theta_rad_vec(kk);
              CIRC_obj       = bh_4bar_kin_circ_CLS( L_1, L_2, L_3, L_4, theta_rad );      
              tf_is_valid_col(kk) =  CIRC_obj.is_valid_pose();
          end
          % done
          tf_is_valid_pose_col = tf_is_valid_col;
        end
        %------------------------------------------------------------------   
        function [z_pair_mat, theta_deg_col, tf_is_real] = calc_z_values(obj)
           theta_deg_col = obj.tgt_theta_deg_col;           
           N             = length(theta_deg_col); 
           z_pair_mat    = zeros(N,2);
           tf_is_real    = false(N,1);
           
           L_1 = obj.L1;
           L_2 = obj.L2;
           L_3 = obj.L3;
           L_4 = obj.L4;

           for kk=1:length(theta_deg_col)
               theta_deg = theta_deg_col(kk);
               theta_rad = deg2rad(theta_deg);

               [the_z, the_z_NUM, the_z_DEN] = bh_z_calc_util(L_1, L_2, L_3, L_4, theta_rad, "CALC_Z");  

               z_pair_mat(kk, 1) = the_z(1);
               z_pair_mat(kk, 2) = the_z(2);
               
               tf_is_real(kk) = isreal(the_z);
               if(tf_is_real(kk))
                   z_pair_mat(kk,:) = real(z_pair_mat(kk,:));
               end
           end
        end
        %------------------------------------------------------------------       
        function [phi_pair_rad, alpha_pair_rad] = calc_phi_alpha_rad(obj, theta_rad)
           L_1 = obj.L1;
           L_2 = obj.L2;
           L_3 = obj.L3;
           L_4 = obj.L4;
            
           CIRC_obj       = bh_4bar_kin_circ_CLS( L_1, L_2, L_3, L_4, theta_rad );
           
           [phi_pair_rad, alpha_pair_rad] = CIRC_obj.calc_phi_alpha_rad();
           
           phi_pair_rad   = phi_pair_rad(:);   % a COLUMN
           alpha_pair_rad = alpha_pair_rad(:); % a COLUMN
        end
        %------------------------------------------------------------------
        function [z_pair, num_pair, den_pair] = calc_z(obj, theta_rad)
           L_1 = obj.L1;
           L_2 = obj.L2;
           L_3 = obj.L3;
           L_4 = obj.L4;
           
           % INTENTION:
           %   Z will be complex IFFF we have an INvalid pose.
           [z_pair, num_pair, den_pair] = bh_z_calc_util(L_1, L_2, L_3, L_4, theta_rad, "CALC_Z");
        end
        %------------------------------------------------------------------        
        function [phi_rad_value, alpha_rad_value] = get_phi_alpha(obj, theta_rad, act_str)
           % ALLOWED USAGE:
           %  >> [phi_rad, alpha_rad] = get_phi_alpha(obj, theta_rad)
           %  >> [phi_rad, alpha_rad] = get_phi_alpha(obj, theta_rad, "CALC_PHI_AND_ALPHA_VAL_AND_PHI_IS_POS")
           %---------------------------------------------------------------
           % ATTENTION:
           %  The default usage 
           %   >> [phi_rad, alpha_rad] = get_phi_alpha(obj, theta_rad)
           %  will return the machine pose angles that correspond to PHI>0.
           %  IFFF it's possible that the machine has 2 positive PHI poses
           %  then PHI_1  (ie the solution from Z_1) will be used.
           %---------------------------------------------------------------
           L_1 = obj.L1;
           L_2 = obj.L2;
           L_3 = obj.L3;
           L_4 = obj.L4;
        
           [phi_pair_rad, alpha_pair_rad] = calc_phi_alpha_rad(obj, theta_rad);
           
           if(2==nargin)
                if( phi_pair_rad(1) >=0 )
                  alpha_rad_value = alpha_pair_rad(1);
                  phi_rad_value   =   phi_pair_rad(1);
                else
                  alpha_rad_value = alpha_pair_rad(2);
                  phi_rad_value   =   phi_pair_rad(2);
                end
           else           
                  alpha_rad_value = alpha_pair_rad(1);
                  phi_rad_value   =   phi_pair_rad(1);
           end
        end
        %------------------------------------------------------------------
        function plot_valid_theta(obj, hax_vec)
            % ALLOWED USAGE:
            %  >> plot_valid_theta(obj)
            %  >> plot_valid_theta(obj, hax_vec)
           if(1==nargin)
               figure;
               hax_vec(1) = subplot(1,2,1);  % round
               hax_vec(2) = subplot(1,2,2);  % tf flat
           end
            
           N = length(obj.tgt_theta_deg_col);
           c_rgb = jet(N);

           x     = cosd(obj.tgt_theta_deg_col);
           y     = sind(obj.tgt_theta_deg_col);
           x(~obj.tf_tgt_is_valid_pose_col) = NaN;
           y(~obj.tf_tgt_is_valid_pose_col) = NaN;
               
           scatter(hax_vec(1),x, y, [], c_rgb);
           xlim(hax_vec(1),[-1,1]);
           ylim(hax_vec(1),[-1,1]);
           
           axis(hax_vec(1),'equal');
           grid(hax_vec(1),'on');
           title(hax_vec(1), 'Valid angular values for \theta');
           
           th_deg = obj.tgt_theta_deg_col;
           tf_vec = double( obj.tf_tgt_is_valid_pose_col );
           
           %th_deg(~obj.tf_tgt_is_valid_pose_col) = NaN;
           %tf_vec(~obj.tf_tgt_is_valid_pose_col)  = NaN;
           
           scatter(hax_vec(2),th_deg, tf_vec, [], c_rgb);
              hold(hax_vec(2), 'on');
              plot(hax_vec(2), th_deg, tf_vec, '-k');
              axis(hax_vec(2),'tight')   
             %xlim(hax_vec(2), [min(obj.tgt_theta_deg_col), max(obj.tgt_theta_deg_col)]);
           ylim(hax_vec(2), [-0.25, 1.25])
           grid(hax_vec(2),'on');
           title(hax_vec(2), 'Valid angular values for \theta');
           xlabel(hax_vec(2), '\theta (deg)')
        end
        %------------------------------------------------------------------                        
        function plot(obj, hax)
            % ALLOWED USAGE:
            %  >> plot(obj)
            %  >> plot(obj, hax_vec)
           if(1==nargin)
               figure;
               hax(1) = subplot(2,1,1);
               hax(2) = subplot(2,1,2);
           end
           
           theta_rad = obj.theta_rad;
          
           [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, "FIRST");      
           LOC_plot(hax(1),A,B1,B2,C,D);
           title(hax(1),"the Z1 solution");
           
           [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, "SECOND");      
           LOC_plot(hax(2),A,B1,B2,C,D);
           title(hax(2),"the Z2 solution");
        end
        %------------------------------------------------------------------
        function plot_is_valid(obj, hax, hax_round)
            % ALLOWED USAGE:
            %  >> plot_is_valid(obj)
            %  >> plot_is_valid(obj, hax)
           if(1==nargin)
               figure;
               hax = axes;
           end
           N      = length(obj.tgt_theta_deg_col);
           c_rgb  = jet(N);

           x      = obj.tgt_theta_deg_col;
           y      = double(obj.tf_tgt_is_valid_pose_col);
   
           scatter(hax, x, y, [], c_rgb); 
           hold(hax, 'on');
           plot(hax, x, y, '-k');
           grid(hax, 'on');
           axis(hax,'tight');   
           ylim(hax, [-0.25, 1.25])
           xlabel(hax,'\theta (degs)');
           ylabel(hax,'Is VALID pose');
           title(hax,'Where does the machine have a VALID pose ?')
        end
        %------------------------------------------------------------------
        function plot_phi_pos(obj, hax)
            % ALLOWED USAGE:
            %  >> plot_phi_pos(obj)
            %  >> plot_phi_pos(obj, hax_vec)
           if(1==nargin)
               figure;
               hax(1) = axes;
           end
           
           theta_rad = obj.theta_rad;
          
           [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, "PHI_POSITIVE_POSE");      
           LOC_plot(hax,A,B1,B2,C,D);
           title(hax(1),"the \phi positive solution");
           
        end
        %------------------------------------------------------------------
        function plot3(obj)
           if(1==nargin)
               figure;
               hax(1) = subplot(2,1,1);
               hax(2) = subplot(2,1,2);
           end
           
           my_stuff_T = LOC_do_theta_deg_sweep(obj);
           
           theta_deg_vec = my_stuff_T.theta_deg_vec;
           phi_1_deg_vec = my_stuff_T.phi_1_deg_vec;
           phi_2_deg_vec = my_stuff_T.phi_2_deg_vec;
           alp_1_deg_vec = my_stuff_T.alp_1_deg_vec;
           alp_2_deg_vec = my_stuff_T.alp_2_deg_vec;
           
           LOC_plot3(hax(1),theta_deg_vec, phi_1_deg_vec, alp_1_deg_vec);
           LOC_plot3(hax(2),theta_deg_vec, phi_2_deg_vec, alp_2_deg_vec);
        end
        %------------------------------------------------------------------     
        function plot_angles_flat(obj, hax)
            % ALLOWED USAGE:
            %  >> plot_angles_flat(obj)
            %  >> plot_angles_flat(obj, hax_vec)
           if(1==nargin)
               figure;
               hax(1) = subplot(2,3,1);
               hax(2) = subplot(2,3,2);
               hax(3) = subplot(2,3,3);
               hax(4) = subplot(2,3,4);
               hax(5) = subplot(2,3,5);
               hax(6) = subplot(2,3,6);             
           end
           
           my_stuff_T = LOC_do_theta_deg_sweep(obj);
           
           theta_deg_vec = my_stuff_T.theta_deg_vec;
           phi_1_deg_vec = my_stuff_T.phi_1_deg_vec;
           phi_2_deg_vec = my_stuff_T.phi_2_deg_vec;
           alp_1_deg_vec = my_stuff_T.alp_1_deg_vec;
           alp_2_deg_vec = my_stuff_T.alp_2_deg_vec;
           
           LOC_plot_angles(hax(1:3),theta_deg_vec, phi_1_deg_vec, alp_1_deg_vec);
           LOC_plot_angles(hax(4:6),theta_deg_vec, phi_2_deg_vec, alp_2_deg_vec);
        end
        %------------------------------------------------------------------ 
        function plot_flat_det_phi_alpha(obj, hax)
            % ALLOWED USAGE:
            %  >> plot_flat_det_phi_alpha(obj)
            %  >> plot_flat_det_phi_alpha(obj, hax_vec)
           if(1==nargin)
               figure;
               hax(1) = subplot(2,1,1);
               hax(2) = subplot(2,1,2);
           end
           
           my_stuff_T = LOC_do_theta_deg_sweep(obj);
           
           theta_deg_vec = my_stuff_T.theta_deg_vec;
           phi_1_deg_vec = my_stuff_T.phi_1_deg_vec;
           phi_2_deg_vec = my_stuff_T.phi_2_deg_vec;
           alp_1_deg_vec = my_stuff_T.alp_1_deg_vec;
           alp_2_deg_vec = my_stuff_T.alp_2_deg_vec;
            
           N = length(theta_deg_vec);

           c_rgb = jet(N);

           AX = hax(1);
           scatter(AX,phi_1_deg_vec, alp_1_deg_vec, [], c_rgb)
               grid(AX, 'on');
               %axis(AX,'equal')
               axis(AX,'tight');
               xlabel(AX,"\phi (deg)");
               ylabel(AX,"\alpha (deg)");
               title(AX, "\phi vs \alpha SOLUTION #1")
           AX = hax(2);
           scatter(AX,phi_2_deg_vec, alp_2_deg_vec, [], c_rgb)
               grid(AX, 'on');
               %axis(AX,'equal')
               axis(AX,'tight');
               xlabel(AX,"\phi (deg)");
               ylabel(AX,"\alpha (deg)");
               title(AX, "\phi vs \alpha SOLUTION #2")
    
          my_det = @(phi, alpha) sind(phi).*cosd(alpha) - cosd(phi).*sind(alpha);   
          
          AX = hax(1);
          hold(AX,'on')   
          %my_lim_vec = [min(phi_1_deg_vec), max(phi_1_deg_vec), min(alp_1_deg_vec), max(alp_1_deg_vec)];
          my_lim_vec = [-180, 180,  -180, 180];
          fcontour(AX, my_det, my_lim_vec, "LevelList",[0],"LineWidth",3)
          
          AX = hax(2);
          hold(AX,'on')   
          %my_lim_vec = [min(phi_2_deg_vec), max(phi_2_deg_vec), min(alp_2_deg_vec), max(alp_2_deg_vec)];
          my_lim_vec = [-180, 180,  -180, 180];
          fcontour(AX, my_det, my_lim_vec, "LevelList",[0],"LineWidth",3)

        end
        %------------------------------------------------------------------         
        function plot_angles_round(obj, hax)
        %  >> plot_angles_round(obj)
        %  >> plot_angles_round(obj, hax_vec)
        %  >> plot_angles_round(obj, {NaN, NaN, htheta2, NaN, NaN, NaN})
        %--------------------------------------------------------------
        % the REQUIRED order of any axes handles MUST be:
        %            1        2          3          4        5          6      
        %   [h_theta_1, h_phi_1, h_alpha_1, h_theta_2, h_phi_2, h_alpha_2]
            
           if(1==nargin)
               figure;
               hax(1) = subplot(2,3,1);
               hax(2) = subplot(2,3,2);
               hax(3) = subplot(2,3,3);
               hax(4) = subplot(2,3,4);
               hax(5) = subplot(2,3,5);
               hax(6) = subplot(2,3,6);    
           else
               assert(6==length(hax),'###_ATTENTION:  your hax vector MUST have 6 elements')
           end
           
           my_stuff_T = LOC_do_theta_deg_sweep(obj);
           N          = length(obj.tgt_theta_deg_col);
           c_rgb      = jet(N);
           
           % so we have a vector of axes handles
           if( ~iscell(hax) )
               LOC_plot_angle_round(hax(1), c_rgb, my_stuff_T.theta_deg_vec);
               LOC_plot_angle_round(hax(2), c_rgb, my_stuff_T.phi_1_deg_vec);
               LOC_plot_angle_round(hax(3), c_rgb, my_stuff_T.alp_1_deg_vec);
               LOC_plot_angle_round(hax(4), c_rgb, my_stuff_T.theta_deg_vec);
               LOC_plot_angle_round(hax(5), c_rgb, my_stuff_T.phi_2_deg_vec);
               LOC_plot_angle_round(hax(6), c_rgb, my_stuff_T.alp_2_deg_vec);

               title(hax(1), "\theta");
               title(hax(2), "\phi_1");
               title(hax(3), "\alpha_1");
               title(hax(4), "\theta");
               title(hax(5), "\phi_2");           
               title(hax(6), "\alpha_2");
           else
                % hax is a cell array 
                if( isgraphics(hax{1}) )
                   LOC_plot_angle_round(hax{1}, c_rgb, my_stuff_T.theta_deg_vec);
                   title(hax{1}, "\theta");
                end
                
                if( isgraphics(hax{2}) )
                  LOC_plot_angle_round(hax{2}, c_rgb, my_stuff_T.phi_1_deg_vec);
                  title(hax{2}, "\phi_1");
                end
                
                if( isgraphics(hax{3}) )
                  LOC_plot_angle_round(hax{3}, c_rgb, my_stuff_T.alp_1_deg_vec);
                  title(hax{3}, "\alpha_1");
                end
                
                if( isgraphics(hax{4}) )
                  LOC_plot_angle_round(hax{4}, c_rgb, my_stuff_T.theta_deg_vec);
                  title(hax{4}, "\theta");
                end
                
                if( isgraphics(hax{5}) )
                  LOC_plot_angle_round(hax{5}, c_rgb, my_stuff_T.phi_2_deg_vec);
                  title(hax{5}, "\phi_2");      
                end
                
                if( isgraphics(hax{6}) )
                  LOC_plot_angle_round(hax{6}, c_rgb, my_stuff_T.alp_2_deg_vec);
                  title(hax{6}, "\alpha_2");
                end
           end % ~iscell()
        end %function plot_angles_round()
        %------------------------------------------------------------------ 
        function plot_z(obj, hax)
            %  >> plot_z(obj)
            %  >> plot_z(obj, hax_vec)
           if(1==nargin)
               figure;
               hax = axes;
           end
           
           my_stuff_T = LOC_do_theta_deg_sweep(obj);
           
%            N = length(obj.tgt_theta_deg_col);
%            c_rgb = jet(N);
           
           X = my_stuff_T.theta_deg_vec;
          Z1 = my_stuff_T.z1_vec;
          Z2 = my_stuff_T.z2_vec;           
           
          sgnLog_Z1 = sign(Z1) .* log10(abs(Z1));
          sgnLog_Z2 = sign(Z2) .* log10(abs(Z2));
          
          plot(hax,X, sgnLog_Z1, '-r.', 'LineWidth',2);
             hold(hax,'on');
          plot(hax,X, sgnLog_Z2, '-b.', 'LineWidth',2);
             grid(hax, 'on');
             axis(hax,'tight');
             xlabel(hax,"\theta (deg)");
             ylabel(hax,"signed Log10( abs(Z))");
             legend(hax,["signed Log10( abs(Z1))","signed Log10( abs(Z2))"])
             title(hax,'signed Log10( abs(Z))')
        end
        %------------------------------------------------------------------ 
        function plot_z_num_den(obj, hax_vec)
            %  >> plot_z_num_den(obj)
            %  >> plot_z_num_den(obj, hax_vec)
           if(1==nargin)
               figure;
               hax_vec(1) = subplot(2,1,1);
               hax_vec(2) = subplot(2,1,2);
           end
        
          % clear the axes 
%           cla(hax_vec(1));
%           cla(hax_vec(2));
          
          % compute the Z1 and Z2 things
          theta_deg_col = obj.tgt_theta_deg_col;
          theta_rad_col = deg2rad( theta_deg_col );
          N             = length(theta_deg_col); 
          the_num_mat = zeros(2,N);
          the_den_mat = zeros(2,N);
          
          for kk=1:N
              [~, the_num, the_den] = bh_z_calc_util(obj.L1, ...
                                                     obj.L2, ...
                                                     obj.L3, ...
                                                     obj.L4, ...
                                                     theta_rad_col(kk), ...
                                                     "CALC_Z");
              if(~isreal(the_num(1)))                                   
                  the_num(1) = NaN;
              end
              if(~isreal(the_num(2)))                                   
                  the_num(2) = NaN;
              end
              if(~isreal(the_den(1)))                                   
                  the_den(1) = NaN;
              end
              if(~isreal(the_den(2)))                                   
                  the_den(2) = NaN;
              end
              
              the_num_mat(:,kk) =  the_num(:);
              the_den_mat(:,kk) =  the_den(:);
          end
         % plot the Z1 stuff
         THE_N  = the_num_mat(1,:);
         THE_D  = the_den_mat(1,:);
         THE_AX = hax_vec(1);
         
         plot(THE_AX, theta_deg_col,  THE_N, '-r.');
           hold(THE_AX,'on');
         plot(THE_AX, theta_deg_col,  THE_D, '-b.');  
           grid(THE_AX,'on');
           xlabel(THE_AX,'theta (deg)'); 
         axis(THE_AX,'tight');
         
         yline(THE_AX,0,'-g', 'LineWidth',2); 
         legend(THE_AX,["The NUMER","The DENOM","zero"],'Location', 'best');
         title(THE_AX, 'Z1, Numerator and Denominator') 
         
         % plot the Z2 stuff
         THE_N  = the_num_mat(2,:);
         THE_D  = the_den_mat(2,:);
         THE_AX = hax_vec(2);
         
         plot(THE_AX, theta_deg_col,  THE_N, '-r.');
           hold(THE_AX,'on');
         plot(THE_AX, theta_deg_col,  THE_D, '-b.');  
           grid(THE_AX,'on');
           xlabel(THE_AX,'theta (deg)'); 
           axis(THE_AX,'tight');
           
           yline(THE_AX,0,'-g', 'LineWidth',2); 
           legend(THE_AX,["The NUMER","The DENOM","zero"],'Location', 'best');
         title(THE_AX, 'Z2, Numerator and Denominator') 
         
        end 
       %------------------------------------------------------------------ 
       function plot_circles(obj, theta_deg, hax)
            %  >> plot_circles(obj, theta_deg)
            %  >> plot_circles(obj, theta_deg, hax_vec)
           if(1==nargin)
               theta_deg = 45;
               figure;
               hax = axes;
           elseif(2==nargin)
               figure;
               hax = axes;
           end
                      
          % clear the axes 
          cla(hax); 
          
          % set the HOLD state
          hold(hax,'on') 
           
          % Create the L1 and L3 polygons circle
          N_POINTS = 720;
          L1_circ = nsidedpoly(N_POINTS,'Radius',obj.L1, 'Center', [0,0]);
          L3_circ = nsidedpoly(N_POINTS,'Radius',obj.L3, 'Center', [obj.L4,0]);
          
          % plot the L1 and L3 polygons
          L1_circ.plot('Parent', hax, 'FaceColor', 'red');
          L3_circ.plot('Parent', hax, 'FaceColor', 'blue');
          
          % create the L2 polygon
          xc = obj.L1*cosd(theta_deg);
          yc = obj.L1*sind(theta_deg);
          L2_circ = nsidedpoly(N_POINTS,'Radius',obj.L2, 'Center', [xc,yc]);
          
          % plot the L2 polygon
          L2_circ.plot('Parent', hax, 'FaceColor', 'green');
          
          plot(hax,xc,yc,'k.','MarkerSize',15)
          
          grid(hax,'on')
          
          % set the axes limits
          the_ymax = max([   (obj.L1+obj.L2),    obj.L3]);
          the_ymin = min([-1*(obj.L1+obj.L2), -1*obj.L3]);
          
          the_xmin = min( [-1*( obj.L1 + obj.L2), (obj.L4-obj.L3)]);
          
          the_xmax = max([(obj.L4+obj.L3),  (obj.L1+obj.L2)]);
          
          
          % make the axes scales equal
          axis(hax,'equal');

          xlim(hax,[the_xmin, the_xmax]);
          ylim(hax,[the_ymin, the_ymax]);
          
          % add the TEXT
          text(hax,0,0,     'C','FontSize', 16,'FontWeight', 'Bold');
          text(hax,obj.L4,0,'D','FontSize', 16,'FontWeight', 'Bold');
          
          xA = obj.L1*cosd(theta_deg);
          yA = obj.L1*sind(theta_deg);
          text(hax,xA,yA,'A','FontSize', 16,'FontWeight', 'Bold');
          
       end
        %------------------------------------------------------------------ 
        function animate_simple(obj)
           if(1==nargin)
               figure;
               hax(1) = subplot(1,2,1);
               hax(2) = subplot(1,2,2);
           end
          
           my_stuff_T = LOC_do_theta_deg_sweep(obj);
           
           theta_deg_vec = my_stuff_T.theta_deg_vec;
           phi_1_deg_vec = my_stuff_T.phi_1_deg_vec;
           phi_2_deg_vec = my_stuff_T.phi_2_deg_vec;
           alp_1_deg_vec = my_stuff_T.alp_1_deg_vec;
           alp_2_deg_vec = my_stuff_T.alp_2_deg_vec;
           
           THE_YMAX =    max([obj.L3, obj.L1]);
           THE_YMIN = -1*THE_YMAX;
           
           THE_XMIN = -1.2*obj.L1;
           THE_XMAX =  obj.L4 + 1.2*max( [obj.L3*cosd(phi_1_deg_vec); ...
                                          obj.L3*cosd(phi_2_deg_vec)] );
                      
           for kk=1:length(obj.tgt_theta_deg_col)
               if(~obj.tf_tgt_is_valid_pose_col(kk))
                   continue
               end
               
               if(~ishghandle(hax(1)) | ~ishghandle(hax(1)) )
                   break
               end
               
               theta_deg = theta_deg_vec(kk);
               theta_rad = deg2rad(theta_deg);

               [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, "FIRST");      
               LOC_plot_for_animate(hax(1),A,B1,B2,C,D);
                              
               [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, "SECOND");      
               LOC_plot_for_animate(hax(2),A,B1,B2,C,D);
               
               for aa=1:2
                  grid(hax(aa), 'on');
                  axis(hax(aa),'equal');
                  xlim(hax(aa), [THE_XMIN, THE_XMAX]);
                  ylim(hax(aa), [THE_YMIN, THE_YMAX]);
                  %    xlabel(hax,"X");
                  %    ylabel(hax,"Y")
               end
               
               drawnow
           end
        end
        %------------------------------------------------------------------     
        function animate_for_app(obj, z_sol_str, app)
            % ALLOWED USAGE:
            %  >> animate_for_app(obj, "FIRST",  app_OBJ)
            %  >> animate_for_app(obj, "SECOND", app_OBJ)
            
            my_stuff_T = LOC_do_theta_deg_sweep(obj);
            
            z_sol_str = upper(z_sol_str);
            switch(z_sol_str)
                case "FIRST"
                     theta_deg_vec = my_stuff_T.theta_deg_vec;
                     phi_deg_vec   = my_stuff_T.phi_1_deg_vec;
                     alp_deg_vec   = my_stuff_T.alp_1_deg_vec;
                     
                     hax_AN    = app.UIAxes_AN_Z1_AN;
                     hax_theta = app.UIAxes_AN_Z1_THETA;
                     hax_phi   = app.UIAxes_AN_Z1_PHI;
                     hax_alpha = app.UIAxes_AN_Z1_ALPHA;
                     hSL       = app.Slider_AN_Z1;
                     hBUT      = app.But_AN_Z1;
                     hSW       = app.STOP_AN_Z1;
                case "SECOND"
                     theta_deg_vec = my_stuff_T.theta_deg_vec;
                     phi_deg_vec   = my_stuff_T.phi_2_deg_vec;
                     alp_deg_vec   = my_stuff_T.alp_2_deg_vec;
                     
                     hax_AN    = app.UIAxes_AN_Z2_AN;
                     hax_theta = app.UIAxes_AN_Z2_THETA;
                     hax_phi   = app.UIAxes_AN_Z2_PHI;
                     hax_alpha = app.UIAxes_AN_Z2_ALPHA;
                     hSL       = app.Slider_AN_Z2;
                     hBUT      = app.But_AN_Z2;
                     hSW       = app.STOP_AN_Z2;
                otherwise
                    error('###_ERROR:  UNknown Z_SOL_STR !');
            end %switch
                        
            % set the "STOP" animation switch to NO
            hSW.Value              = 'NO';
            
            % what are the user's PAUSE times in milliseconds
            pause_time_for_valid   = app.PAUSE_ValidanglesEditField.Value;
            pause_time_for_invalid = app.PAUSE_INvalidanglesEditField.Value;
            
            % what is the display frame frequency
            N_display_frame        = str2num(app.DISPLAY_EVERY_N_FRAMES_Knob.Value);
                       
            % compute display LIMITS
            [the_xlim, the_ylim] = LOC_compute_XY_display_limits_for_animate(obj, my_stuff_T, z_sol_str);
            THE_XMIN = the_xlim(1);
            THE_XMAX = the_xlim(2);
            THE_YMIN = the_ylim(1);
            THE_YMAX = the_ylim(2);
                                   
            grid(hax_AN, 'on');
            axis(hax_AN,'equal');
            xlim(hax_AN, [THE_XMIN, THE_XMAX]);
            ylim(hax_AN, [THE_YMIN, THE_YMAX]);
 
            % disable the ANIMATE pushbutton
            hBUT.Enable = 'off';
            hBUT.BackgroundColor = [0.96,0.96,0.96];
            hSL.Enable = 'off';
            
            % configure the HOLD attribute of our target AXES
            hold(hax_theta, 'on');
            hold(hax_phi,   'on');
            hold(hax_alpha, 'on');
            hold(hax_AN, 'on');
            
            % setup some ANIMATEDLINES for the 3 ANGULAR axes
            hL_theta = findobj(hax_theta, 'TAG', 'TAG_LINE_THETA');
            hL_phi   = findobj(hax_phi  , 'TAG', 'TAG_LINE_PHI');
            hL_alpha = findobj(hax_alpha, 'TAG', 'TAG_LINE_ALPHA');

            delete(hL_theta);
            delete(hL_phi);
            delete(hL_alpha);  

            hL_theta = animatedline(hax_theta, 'Color','r','LineWidth',3,'TAG', 'TAG_LINE_THETA');
            hL_phi   = animatedline(hax_phi  , 'Color','b','LineWidth',3,'TAG', 'TAG_LINE_PHI');
            hL_alpha = animatedline(hax_alpha, 'Color','g','LineWidth',3,'TAG', 'TAG_LINE_ALPHA');
            
            % setup some ANIMATEDLINES for the 4BAR axes
            hL_4bar_theta = findobj(hax_AN, 'TAG', 'TAG_4BAR_LINE_THETA');
            hL_4bar_phi   = findobj(hax_AN, 'TAG', 'TAG_4BAR_LINE_PHI');
            hL_4bar_alpha = findobj(hax_AN, 'TAG', 'TAG_4BAR_LINE_ALPHA');
                  
            delete(hL_4bar_theta);
            delete(hL_4bar_phi);
            delete(hL_4bar_alpha);  
                  
            hL_4bar_theta = animatedline(hax_AN, 'Color','r','LineWidth',3,'TAG', 'TAG_4BAR_LINE_THETA');
            hL_4bar_phi   = animatedline(hax_AN, 'Color','b','LineWidth',3,'TAG', 'TAG_4BAR_LINE_PHI');
            hL_4bar_alpha = animatedline(hax_AN, 'Color','g','LineWidth',3,'TAG', 'TAG_4BAR_LINE_ALPHA');
            
           a = tic; 
           % OK - cycle through the target theta values 
           for kk=1:length(obj.tgt_theta_deg_col)
               
               % should I STOP the animation ?
               if( "YES" == string(hSW.Value) )
                   break
               end
               
               % should I update the display ?
               if(0 ~= mod(kk,N_display_frame)  &&  kk~=1)
                   continue
               end     
               
               % what's THETA               
               tgt_theta_deg = obj.tgt_theta_deg_col(kk);
               
               hSL.Value     = tgt_theta_deg;
               
               x             = cosd(tgt_theta_deg);  
               y             = sind(tgt_theta_deg);

               % update our THETA axes plot
               clearpoints(hL_theta);
               addpoints(hL_theta,[0,x],[0,y])
                              
               if(~obj.tf_tgt_is_valid_pose_col(kk))
                   clearpoints(hL_phi);
                   clearpoints(hL_alpha);
                   set(hax_AN,'Color',[1,0.8,0.5]);
                   drawnow
                   pause(pause_time_for_invalid/1000);
                   continue
               else
                   set(hax_AN,'Color',[1,1,1]);
                   phi_deg   =   phi_deg_vec(kk); 
                   alpha_deg =   alp_deg_vec(kk);
                   theta_deg = theta_deg_vec(kk); 
                   theta_rad = deg2rad(theta_deg);  
                                 
                   [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, z_sol_str);
                   
                   %update the 4bar axes
                   clearpoints(hL_4bar_theta);
                   clearpoints(hL_4bar_phi);
                   clearpoints(hL_4bar_alpha);
                  
                   addpoints(hL_4bar_theta, [ C(1),  A(1)], [ C(2),  A(2)]); 
                   addpoints(hL_4bar_alpha, [ A(1), B1(1)], [ A(2), B1(2)]); 
                   addpoints(hL_4bar_phi,   [B2(1),  D(1)], [B2(2),  D(2)]); 
                   
                   % update our PHI axes plot
                   x = cosd(phi_deg);  y = sind(phi_deg);
                   clearpoints(hL_phi);
                   addpoints(hL_phi,[0,x],[0,y]);

                   % update our ALPHA axes plot
                   x = cosd(alpha_deg);  y = sind(alpha_deg);
                   clearpoints(hL_alpha);
                   addpoints(hL_alpha,[0,x],[0,y]);
                   
%                    if(kk==1 | 0==mod(kk,3))
%                        drawnow
%                    end
                   
                  if(0==mod(kk,N_display_frame))
                   drawnow
                  end
                  
                  if(pause_time_for_valid > 0)
                      pause(pause_time_for_valid/1000);
                  end
               end
               %fprintf('\n ... STEP %3d of %d ',kk,length(obj.tgt_theta_deg_col));               
           end % for kk=1:length(obj.tgt_theta_deg_col)
           
           % RE enablke some of the controls
            hBUT.Enable = 'on';
            hBUT.BackgroundColor = [1,1,0.07];
            
            hSL.Enable = 'on';

        end
        %------------------------------------------------------------------
        function slider_for_app(obj, z_sol_str, app)
            % ALLOWED USAGE:
            %  >> slider_for_app(obj, "FIRST",  app_OBJ)
            %  >> slider_for_app(obj, "SECOND", app_OBJ)
               
            my_stuff_T = LOC_do_theta_deg_sweep(obj);

            
            z_sol_str = upper(z_sol_str);
            switch(z_sol_str)
                case "FIRST"    
                     phi_deg_vec   = my_stuff_T.phi_1_deg_vec;

                     hax_AN    = app.UIAxes_AN_Z1_AN;
                     hax_theta = app.UIAxes_AN_Z1_THETA;
                     hax_phi   = app.UIAxes_AN_Z1_PHI;
                     hax_alpha = app.UIAxes_AN_Z1_ALPHA;
                     hSL       = app.Slider_AN_Z1;
                     hBUT      = app.But_AN_Z1;
                     hEd       = app.theta_edit_Z1;
                     htxt.phi  = app.Text_phi_Z1;
                     htxt.alpha= app.Text_alpha_Z1;
                     htxt.theta= app.Text_theta_Z1;
                case "SECOND"                     
                     phi_deg_vec   = my_stuff_T.phi_2_deg_vec;

                     hax_AN    = app.UIAxes_AN_Z2_AN;
                     hax_theta = app.UIAxes_AN_Z2_THETA;
                     hax_phi   = app.UIAxes_AN_Z2_PHI;
                     hax_alpha = app.UIAxes_AN_Z2_ALPHA;
                     hSL       = app.Slider_AN_Z2;
                     hBUT      = app.But_AN_Z2;
                     hEd       = app.theta_edit_Z2;
                     htxt.phi  = app.Text_phi_Z2;
                     htxt.alpha= app.Text_alpha_Z2;
                     htxt.theta= app.Text_theta_Z2;
                otherwise
                    error('###_ERROR:  UNknown Z_SOL_STR !');
            end %switch
            
            [the_xlim, the_ylim] = LOC_compute_XY_display_limits_for_animate(obj, my_stuff_T, z_sol_str);
            THE_XMIN = the_xlim(1);
            THE_XMAX = the_xlim(2);
            THE_YMIN = the_ylim(1);
            THE_YMAX = the_ylim(2);
                
            grid(hax_AN, 'on');
            axis(hax_AN,'equal');
            xlim(hax_AN, [THE_XMIN, THE_XMAX]);
            ylim(hax_AN, [THE_YMIN, THE_YMAX]);
 
            theta_deg = hSL.Value; 
            theta_rad = deg2rad(theta_deg);

            % configure the HOLD attribute of our target AXES
            hold(hax_theta, 'on');
            hold(hax_phi,   'on');
            hold(hax_alpha, 'on');
            hold(hax_AN, 'on');
            
            % setup some ANIMATEDLINES for the 3 ANGULAR axes
            hL_theta = findobj(hax_theta, 'TAG', 'TAG_LINE_THETA');
            hL_phi   = findobj(hax_phi  , 'TAG', 'TAG_LINE_PHI');
            hL_alpha = findobj(hax_alpha, 'TAG', 'TAG_LINE_ALPHA');

            delete(hL_theta);
            delete(hL_phi);
            delete(hL_alpha);  

            hL_theta = animatedline(hax_theta, 'Color','r','LineWidth',3,'TAG', 'TAG_LINE_THETA');
            hL_phi   = animatedline(hax_phi  , 'Color','b','LineWidth',3,'TAG', 'TAG_LINE_PHI');
            hL_alpha = animatedline(hax_alpha, 'Color','g','LineWidth',3,'TAG', 'TAG_LINE_ALPHA');
            
            % setup some ANIMATEDLINES for the 4BAR axes
            hL_4bar_theta = findobj(hax_AN, 'TAG', 'TAG_4BAR_LINE_THETA');
            hL_4bar_phi   = findobj(hax_AN, 'TAG', 'TAG_4BAR_LINE_PHI');
            hL_4bar_alpha = findobj(hax_AN, 'TAG', 'TAG_4BAR_LINE_ALPHA');
                  
            delete(hL_4bar_theta);
            delete(hL_4bar_phi);
            delete(hL_4bar_alpha);  
                  
            hL_4bar_theta = animatedline(hax_AN, 'Color','r','LineWidth',3,'TAG', 'TAG_4BAR_LINE_THETA');
            hL_4bar_phi   = animatedline(hax_AN, 'Color','b','LineWidth',3,'TAG', 'TAG_4BAR_LINE_PHI');
            hL_4bar_alpha = animatedline(hax_AN, 'Color','g','LineWidth',3,'TAG', 'TAG_4BAR_LINE_ALPHA');
                      
            % update our THETA axes plot           
            x         = cosd(theta_deg);  
            y         = sind(theta_deg);
            clearpoints(hL_theta);
            addpoints(hL_theta,[0,x],[0,y])
            
            if( obj.is_valid_pose(theta_rad)  )
               set(hax_AN,'Color',[1,1,1]);
               [phi_pair_rad, alpha_pair_rad] = obj.calc_phi_alpha_rad(theta_rad);
               
               if(z_sol_str=="FIRST")
                   phi_deg   =   rad2deg(  phi_pair_rad(1)); 
                   alpha_deg =   rad2deg(alpha_pair_rad(1)); 
               else
                   phi_deg   =   rad2deg(  phi_pair_rad(2)); 
                   alpha_deg =   rad2deg(alpha_pair_rad(2)); 
               end

               [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, z_sol_str);      

               %update the 4bar axes
               clearpoints(hL_4bar_theta);
               clearpoints(hL_4bar_phi);
               clearpoints(hL_4bar_alpha);
                  
               addpoints(hL_4bar_theta, [ C(1),  A(1)], [ C(2),  A(2)]); 
               addpoints(hL_4bar_alpha, [ A(1), B1(1)], [ A(2), B1(2)]); 
               addpoints(hL_4bar_phi,   [B2(1),  D(1)], [B2(2),  D(2)]); 
                   
               % update our PHI axes plot
               x = cosd(phi_deg);  y = sind(phi_deg);
               clearpoints(hL_phi);
               addpoints(hL_phi,[0,x],[0,y]);

               % update our ALPHA axes plot
               x = cosd(alpha_deg);  y = sind(alpha_deg);
               clearpoints(hL_alpha);
               addpoints(hL_alpha,[0,x],[0,y]);

               drawnow 
               
               htxt.phi.Text   = sprintf("%.2f",phi_deg);
               htxt.alpha.Text = sprintf("%.2f",alpha_deg);
               htxt.theta.Text = sprintf("%.2f",theta_deg);
               hEd.Value       = str2num(num2str(theta_deg,"%.2f"));
            else
                set(hax_AN,'Color',[1,0.8,0.5]);
                htxt.phi.Text   = "-------";
                htxt.alpha.Text = "-------";
                htxt.theta.Text = sprintf("%.2f",theta_deg);
                hEd.Value       = str2num(num2str(theta_deg,"%.2f"));
             end
                                         
        end
        %------------------------------------------------------------------     
        function animate_with_circles_for_app(obj, z_sol_str, app)
            % ALLOWED USAGE:
            %  >> animate_with_circles_for_app(obj, "FIRST",  app_OBJ)
            %  >> animate_with_circles_for_app(obj, "SECOND", app_OBJ)
            
            my_stuff_T = LOC_do_theta_deg_sweep(obj);
            
            z_sol_str = upper(z_sol_str);
            switch(z_sol_str)
                case "FIRST"
                     theta_deg_vec = my_stuff_T.theta_deg_vec;
                     phi_deg_vec   = my_stuff_T.phi_1_deg_vec;
                     alp_deg_vec   = my_stuff_T.alp_1_deg_vec;
                     
                     hax_AN    = app.UIAxes_Z1_CIRC;
                     hax_theta = app.UIAxes_THETA_Z1_CIRC;
                     hBUT      = app.Button_AN_Z1_CIRC;
                     hSW       = app.STOP_AN_Z1_CIRC;
                case "SECOND"
                     theta_deg_vec = my_stuff_T.theta_deg_vec;
                     phi_deg_vec   = my_stuff_T.phi_2_deg_vec;
                     alp_deg_vec   = my_stuff_T.alp_2_deg_vec;
                     
                     hax_AN    = app.UIAxes_Z2_CIRC;
                     hax_theta = app.UIAxes_THETA_Z2_CIRC;
                     hBUT      = app.Button_AN_Z2_CIRC;
                     hSW       = app.STOP_AN_Z2_CIRC;
                otherwise
                    error('###_ERROR:  UNknown Z_SOL_STR !');
            end %switch
            
            % what are the user's PAUSE times in milliseconds
            pause_time_for_valid   = app.PAUSE_ValidanglesEditField.Value;
            pause_time_for_invalid = app.PAUSE_INvalidanglesEditField.Value;
            
            % what is the display frame frequency
            N_display_frame        = str2num(app.DISPLAY_EVERY_N_FRAMES_Knob.Value);

            % set the "STOP" animation switch to NO
            hSW.Value              = 'NO';
            
            % disable the ANIMATE pushbutton
            hBUT.Enable = 'off';
            hBUT.BackgroundColor = [0.96,0.96,0.96];
            hSL.Enable = 'off';
            
            % clear the axes
            cla(hax_AN);
            cla(hax_theta);
            
            % plot the THETA axes
            obj.plot_angles_round({hax_theta,NaN,NaN,NaN,NaN,NaN});
            axis(hax_theta,'equal');
            xlim(hax_theta, [-1.2, 1.2]);
            ylim(hax_theta, [-1.2, 1.2]);
            hold(hax_theta, 'on');

            % configure the HOLD attribute of our target AXES
            hold(hax_AN, 'on');

            % plot the circles
            %--------------------------------------------------------------
            % Create the L1 and L3 polygons circle
            L1_circ = nsidedpoly(720,'Radius',obj.L1, 'Center', [0,0]);
            L3_circ = nsidedpoly(720,'Radius',obj.L3, 'Center', [obj.L4,0]);
          
            % plot the L1 and L3 polygons
            L1_circ.plot('Parent', hax_AN, 'FaceColor', 'red');
            L3_circ.plot('Parent', hax_AN, 'FaceColor', 'blue');
            %--------------------------------------------------------------
                                   
            % set the axes limits
            THE_YMAX = max([   (obj.L1+obj.L2),    obj.L3]);
            THE_YMIN = min([-1*(obj.L1+obj.L2), -1*obj.L3]);          
            THE_XMIN = min( [-1*( obj.L1 + obj.L2), (obj.L4-obj.L3)]);          
            THE_XMAX = max([(obj.L4+obj.L3),  (obj.L1+obj.L2)]);
                       
            grid(hax_AN, 'on');
            axis(hax_AN,'equal');
            xlim(hax_AN, [THE_XMIN, THE_XMAX]);
            ylim(hax_AN, [THE_YMIN, THE_YMAX]);
            
            % setup some ANIMATEDLINES for the THETA axes
            hL_theta = animatedline(hax_theta, 'Color','r','LineWidth',3,'TAG', 'TAG_LINE_THETA');
            
            % setup some ANIMATEDLINES for the 4BAR axes                 
            hL_4bar_theta = animatedline(hax_AN, 'Color','r','LineWidth',3,'TAG', 'TAG_4BAR_LINE_THETA');
            hL_4bar_phi   = animatedline(hax_AN, 'Color','b','LineWidth',3,'TAG', 'TAG_4BAR_LINE_PHI');
            hL_4bar_alpha = animatedline(hax_AN, 'Color','g','LineWidth',3,'TAG', 'TAG_4BAR_LINE_ALPHA');
            
           a = tic; 
           % OK - cycle through the target theta values 
           for kk=1:length(obj.tgt_theta_deg_col)
               
               % should I STOP the animation ?
               if( "YES" == string(hSW.Value) )
                   break
               end
               
               % should I update the display ?
               if(0 ~= mod(kk,N_display_frame)  &&  kk~=1)
                   continue
               end               
               % what's THETA               
               tgt_theta_deg = obj.tgt_theta_deg_col(kk);
               
               hSL.Value     = tgt_theta_deg;
               x             = cosd(tgt_theta_deg);  
               y             = sind(tgt_theta_deg);

               % update our THETA axes plot
               clearpoints(hL_theta);
               addpoints(hL_theta,[0,x],[0,y])
               
               % update the animation plot 
               xA             = obj.L1*cosd(tgt_theta_deg);  
               yA             = obj.L1*sind(tgt_theta_deg);
               if(1==kk)
                    hA_txt = text(hax_AN,xA,yA,'A','FontSize', 16,'FontWeight', 'Bold','Tag','TAG_A_TEXT');
               else
                   hA_txt.Position = [xA,yA,0];
               end
               
               % create the L2 polygon
               xc = obj.L1*cosd(tgt_theta_deg);
               yc = obj.L1*sind(tgt_theta_deg);
               L2_circ = nsidedpoly(720,'Radius',obj.L2, 'Center', [xc,yc]);

               % plot the L2 polygon
               tmp_thing = findobj(hax_AN,'Tag', 'TAG_L2_CIRCLE');
               delete(tmp_thing);
               L2_circ.plot('Parent', hax_AN, 'FaceColor', 'green','Tag','TAG_L2_CIRCLE');
               
               if(~obj.tf_tgt_is_valid_pose_col(kk))
                   drawnow
                   continue
               else
                   phi_deg   =   phi_deg_vec(kk); 
                   alpha_deg =   alp_deg_vec(kk);
                   theta_deg = theta_deg_vec(kk); 
                   theta_rad = deg2rad(theta_deg);  
                                 
                   [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, z_sol_str);
                   
                   %update the 4bar axes
                   clearpoints(hL_4bar_theta);
                   clearpoints(hL_4bar_phi);
                   clearpoints(hL_4bar_alpha);
                  
                   addpoints(hL_4bar_theta, [ C(1),  A(1)], [ C(2),  A(2)]); 
                   addpoints(hL_4bar_alpha, [ A(1), B1(1)], [ A(2), B1(2)]); 
                   addpoints(hL_4bar_phi,   [B2(1),  D(1)], [B2(2),  D(2)]); 
                                                        
                  if(0==mod(kk,N_display_frame))
                   drawnow
                  end
                  
                  if(pause_time_for_valid > 0)
                      pause(pause_time_for_valid/1000);
                  end
                  
               end
               %fprintf('\n ... STEP %3d of %d ',kk,length(obj.tgt_theta_deg_col));               
           end % for kk=1:length(obj.tgt_theta_deg_col)
           
           % RE enablke some of the controls
            hBUT.Enable = 'on';
            hBUT.BackgroundColor = [1,1,0.07];
            
            hSL.Enable = 'on';

        end
        %------------------------------------------------------------------         
    end % methods
    %_#########################################################################
    % PROTECTED methods
    %_#########################################################################
    methods(Access = protected)
        function ALG_is_valid_pose(obj)
            
        end
        %------------------------------------------------------------------  
        %------------------------------------------------------------------  
        %------------------------------------------------------------------  
              
    end % methods(Access = protected)
end % classdef
%_#########################################################################
% Support functions
%_#########################################################################
function [the_xlim, the_ylim] = LOC_compute_XY_display_limits_for_animate(obj, my_stuff_T, z_sol_str)
    z_sol_str = upper(z_sol_str);
    switch(z_sol_str)
        case "FIRST"
             theta_deg_vec = my_stuff_T.theta_deg_vec;
             phi_deg_vec   = my_stuff_T.phi_1_deg_vec;
             alp_deg_vec   = my_stuff_T.alp_1_deg_vec;
        case "SECOND"
             theta_deg_vec = my_stuff_T.theta_deg_vec;
             phi_deg_vec   = my_stuff_T.phi_2_deg_vec;
             alp_deg_vec   = my_stuff_T.alp_2_deg_vec;
        otherwise
            error('###_ERROR:  UNknown Z_SOL_STR !');
    end %switch
            

    % compute display LIMITS
    THE_YMAX = max([obj.L3*sind(phi_deg_vec); obj.L1*sind(theta_deg_vec)]);
    THE_YMIN = min([obj.L3*sind(phi_deg_vec); obj.L1*sind(theta_deg_vec)]);

    THE_XMIN = min([0; (obj.L4 + obj.L3*cosd(phi_deg_vec)); obj.L1*cosd(theta_deg_vec)]);
    THE_XMAX = max([obj.L4; (obj.L4 + obj.L3*cosd(phi_deg_vec)); obj.L1*cosd(theta_deg_vec)]);

    the_xlim = [THE_XMIN, THE_XMAX];
    the_ylim = [THE_YMIN, THE_YMAX];
end
%--------------------------------------------------------------------------
function [A,B1,B2,C,D] = LOC_get_points(obj, theta_rad, action_str)

   [phi_pair_rad, alpha_pair_rad] = calc_phi_alpha_rad(obj, theta_rad);

   switch(upper(action_str))
       case "FIRST"
            THE_PHI   = phi_pair_rad(1);
            THE_ALPHA = alpha_pair_rad(1);
       case "SECOND"
            THE_PHI   = phi_pair_rad(2);
            THE_ALPHA = alpha_pair_rad(2);
       case "PHI_POSITIVE_POSE"
              if(phi_pair_rad(1) >= 0 )
                  THE_PHI   =   phi_pair_rad(1);
                  THE_ALPHA = alpha_pair_rad(1);
              else
                  THE_PHI   =   phi_pair_rad(2);
                  THE_ALPHA = alpha_pair_rad(2);
              end
       otherwise
          error("ERROR:  unknown action string !");     
   end
   
   THE_THETA = theta_rad;

   C  = [0,0];
   A  = [obj.L1*cos(THE_THETA), obj.L1*sin(THE_THETA)];
   B1 = [ (A(1) + obj.L2*cos(THE_ALPHA)), ...
          (A(2) + obj.L2*sin(THE_ALPHA))  ];
   D  = [obj.L4,0];
   B2 = [( D(1) + obj.L3*cos(THE_PHI)), (D(2)+ obj.L3*sin(THE_PHI))];
end
%--------------------------------------------------------------------------
function LOC_plot(hax,A,B1,B2,C,D)

   plot(hax, [C(1), A(1)], [C(2), A(2)], '-r', 'LineWidth',3); 
   
      hold(hax, 'on')
   
   plot(hax, [A(1), B1(1)], [A(2), B1(2)], '-g', 'LineWidth',3); 
   plot(hax, [B2(1), D(1)], [B2(2), D(2)], '-b', 'LineWidth',3); 
   grid(hax, 'on');
   axis(hax,'equal')
   xlabel(hax,"X");
   ylabel(hax,"Y")
   
   hold(hax, 'off')
end
%--------------------------------------------------------------------------
function LOC_plot_for_animate(hax,A,B1,B2,C,D)
% ALLOWED USAGE:
%    >>   LOC_plot_for_animate(hax,A,B1,B2,C,D)

   hL_theta = findobj(hax, 'Tag', 'TAG_4BAR_LINE_THETA');
   hL_phi   = findobj(hax, 'Tag', 'TAG_4BAR_LINE_PHI');
   hL_alpha = findobj(hax, 'Tag', 'TAG_4BAR_LINE_ALPHA');

   hL = [hL_theta, hL_phi, hL_alpha];

    if(isempty(hL_theta))
       plot(hax, [C(1), A(1)], [C(2), A(2)], '-r', 'LineWidth',3,'Tag', 'TAG_4BAR_LINE_THETA'); 

          hold(hax, 'on')

       plot(hax, [A(1), B1(1)], [A(2), B1(2)], '-g', 'LineWidth',3,'Tag', 'TAG_4BAR_LINE_ALPHA'); 
       plot(hax, [B2(1), D(1)], [B2(2), D(2)], '-b', 'LineWidth',3,'Tag', 'TAG_4BAR_LINE_PHI');    
       hold(hax, 'off')
    else
        set(hL,'Visible', 'off');
        set(hL_theta, 'XData', [ C(1),  A(1)],  'YData', [ C(2),  A(2)]);
        set(hL_alpha, 'XData', [ A(1), B1(1)], 'YData',  [ A(2), B1(2)]);
        set(hL_phi,   'XData', [B2(1),  D(1)], 'YData',  [B2(2),  D(2)]);
        set(hL,'Visible', 'on');
    end

end
%--------------------------------------------------------------------------
function LOC_plot_for_animate_v2(hax,A,B1,B2,C,D)
% ALLOWED USAGE:
%    >>   LOC_plot_for_animate_v2(hax,A,B1,B2,C,D)
   cla(hax);
       plot(hax, [C(1), A(1)], [C(2), A(2)], '-r', 'LineWidth',3,'Tag', 'TAG_4BAR_LINE_THETA'); 

          hold(hax, 'on')

       plot(hax, [A(1), B1(1)], [A(2), B1(2)], '-g', 'LineWidth',3,'Tag', 'TAG_4BAR_LINE_ALPHA'); 
       plot(hax, [B2(1), D(1)], [B2(2), D(2)], '-b', 'LineWidth',3,'Tag', 'TAG_4BAR_LINE_PHI'); 
end
%--------------------------------------------------------------------------
function LOC_plot3(hax,theta_deg_vec, phi_deg_vec, alpha_deg_vec)

   scatter3(hax, theta_deg_vec, phi_deg_vec, alpha_deg_vec);
   
   grid(hax, 'on');
   %axis(hax,'equal')
   xlabel(hax,"\theta (deg)");
   ylabel(hax,"\phi (deg)")
   zlabel(hax,"\alpha (deg)")
end
%--------------------------------------------------------------------------
function my_stuff_T = LOC_do_theta_deg_sweep(obj)    
   L_1 = obj.L1;
   L_2 = obj.L2;
   L_3 = obj.L3;
   L_4 = obj.L4;

   theta_deg_vec = NaN(size(obj.tgt_theta_deg_col));
   phi_1_deg_vec = NaN(size(obj.tgt_theta_deg_col));
   phi_2_deg_vec = NaN(size(obj.tgt_theta_deg_col));
   alp_1_deg_vec = NaN(size(obj.tgt_theta_deg_col));
   alp_2_deg_vec = NaN(size(obj.tgt_theta_deg_col));
   z1_vec        = NaN(size(obj.tgt_theta_deg_col));
   z2_vec        = NaN(size(obj.tgt_theta_deg_col));

   for kk=1:length(obj.tgt_theta_deg_col)
       theta_deg = obj.tgt_theta_deg_col(kk);
       theta_rad = deg2rad(theta_deg);

       CIRC_obj       = bh_4bar_kin_circ_CLS( L_1, L_2, L_3, L_4, theta_rad );
      [phi_pair_rad, alpha_pair_rad] = CIRC_obj.calc_phi_alpha_rad();                          
      [z_pair, ~, ~]                 = calc_z(obj, theta_rad);
      %-----------------------------------------------------------------------
      z1_vec(kk)        = z_pair(1);
      z2_vec(kk)        = z_pair(2);

      phi_pair_deg = rad2deg(phi_pair_rad);
      alp_pair_deg = rad2deg(alpha_pair_rad);

      phi_1_deg_vec(kk) = phi_pair_deg(1);
      phi_2_deg_vec(kk) = phi_pair_deg(2);
      alp_1_deg_vec(kk) = alp_pair_deg(1);
      alp_2_deg_vec(kk) = alp_pair_deg(2);

      if(true==obj.tf_tgt_is_valid_pose_col(kk))       
              theta_deg_vec(kk) = theta_deg;
       end                   
   end % for

   % just test for angle ranges
   fh_test_in_range = @(x) assert( max(x)<=180 && min(x)>=-180);
   fh_test_in_range(theta_deg_vec);
   fh_test_in_range(phi_1_deg_vec);
   fh_test_in_range(phi_2_deg_vec);
   fh_test_in_range(alp_1_deg_vec);
   fh_test_in_range(alp_2_deg_vec);

   % assign function OUTPUTS
   my_stuff_T.tgt_theta_deg_vec = obj.tgt_theta_deg_col;
   my_stuff_T.theta_deg_vec = theta_deg_vec;
   my_stuff_T.phi_1_deg_vec = phi_1_deg_vec;
   my_stuff_T.phi_2_deg_vec = phi_2_deg_vec;
   my_stuff_T.alp_1_deg_vec = alp_1_deg_vec;
   my_stuff_T.alp_2_deg_vec = alp_2_deg_vec;
   my_stuff_T.z1_vec        = z1_vec;
   my_stuff_T.z2_vec        = z2_vec;
end        
%--------------------------------------------------------------------------
function LOC_plot_angles(hax_vec,theta_deg_vec, phi_deg_vec, alp_deg_vec)
   
   N = length(theta_deg_vec);
   
   c_rgb = jet(N);

   AX = hax_vec(1);
   %plot(AX, theta_deg_vec, phi_deg_vec, '-k');
   scatter(AX,theta_deg_vec, phi_deg_vec, [], c_rgb)
       grid(AX, 'on');
       %axis(AX,'equal');
       axis(AX,'tight');
       xlabel(AX,"\theta (deg)");
       ylabel(AX,"\phi (deg)");
       title(AX, "\phi");
   AX = hax_vec(2);
   %plot(AX, theta_deg_vec, alp_deg_vec, '-m');
   scatter(AX,theta_deg_vec, alp_deg_vec, [], c_rgb)
       grid(AX, 'on');
       %axis(AX,'equal')
       axis(AX,'tight');
       xlabel(AX,"\theta (deg)");
       ylabel(AX,"\alpha (deg)");
       title(AX, "\alpha");
   AX = hax_vec(3);
   %plot(AX, phi_deg_vec, alp_deg_vec, '-b');
   scatter(AX,phi_deg_vec, alp_deg_vec, [], c_rgb)
       grid(AX, 'on');
       %axis(AX,'equal')
       axis(AX,'tight');
       xlabel(AX,"\phi (deg)");
       ylabel(AX,"\alpha (deg)");
       title(AX, "\phi vs \alpha")
       
%     my_det = @(phi, alpha) sind(phi).*cosd(alpha) - cosd(phi).*sind(alpha);   
%     hold(AX,'on')   
%     fcontour(AX, my_det, [-180 180  -180 180], "LevelList",[0],"LineWidth",3)
%        xlim(AX, [-180, 180]);
%        ylim(AX, [-180, 180]);
    
    % set the display limits   
%     for kk=1:3
%        xlim(hax_vec(kk), [-180, 180]);
%        ylim(hax_vec(kk), [-180, 180]);
%     end
end
%--------------------------------------------------------------------------
function LOC_plot_angle_round(hax, c_rgb, THE_DEG_ANG_VEC)

           N = length(THE_DEG_ANG_VEC);
           
           R = linspace(1,0.5, N);
           R = reshape(R, size(THE_DEG_ANG_VEC));
           
           x     = R .* cosd(THE_DEG_ANG_VEC);
           y     = R .* sind(THE_DEG_ANG_VEC);
               
           scatter(hax,x, y, [], c_rgb);           
           axis(hax,'equal');
           grid(hax,'on');
           xlim(hax,[-1,1]);
           ylim(hax,[-1,1]);

end
%--------------------------------------------------------------------------
