classdef bh_rot_active_B2G_CLS
% BH_ROT_ACTIVE_B2G is an ACTIVELY rotate a Bvec into a new G location
%==========================================================================
% Assume that a body undergoes 3 successive rotations in the following order:
%
%   a.) R1 occurs 1st about some local body axis
%   b.) R2 occurs 2nd about some local body axis
%   c.) R3 occurs 3rd about some local body axis
%
% First consider the PASSIVE rotation which relates the final G vector to
% it's corresponding components in the B frame:
%
%    vB = R3(c) * R2(b) * R1(a) * vG
%    vB = bRg * vG
%
% Therefore we can now write:
%
%  vG =  R1(-a)  *  R2(-b) *  R3(-c) * vB
%  vG = aR1(+a)  * aR2(+b) * aR3(+c) * vB 
%  vG = agRb * vB
%==========================================================================
% Example USAGE:
%  See the script files:
%     a.) bh_test_rot_active_B2G_example_1.m
%     b.) bh_test_rot_active_B2G_example_2.m
%     c.) bh_test_rot_active_B2G_example_3.m
%--------------------------------------------------------------------------

    properties (SetAccess = protected)
        ang_units ; % categorical {'DEGREES', 'RADIANS', 'SYM'}   
        %---------------------------
        num_rotations = 0;
        %---------------------------
        dir_1st; % categorical {'D1X', 'D1Y', 'D1Z', '---'}
        dir_2nd; % categorical {'D2X', 'D2Y', 'D2Z', '---'}
        dir_3rd; % categorical {'D3X', 'D3Y', 'D3Z', '---'}
        %---------------------------
        ang_1st = 0;
        ang_2nd = 0;
        ang_3rd = 0;
    end
    
    methods
        %------------------------------------------------------------------
        function OBJ = bh_rot_active_B2G_CLS(dir_cell_list, dir_ang_list, units_str)
            
            % set the ANG_UNITS propert first
            OBJ.ang_units = units_str;
            
            dir_cell_list = dir_cell_list(:);
            dir_ang_list  = dir_ang_list(:);
            
            % check dir_cell_list
            assert(length(dir_cell_list)==length(dir_ang_list), 'Why different lengths ?');
            validateattributes(dir_cell_list,{'cell'},{'vector'});
            assert(iscellstr(dir_cell_list),'it needs to be a cell list of strings')
            assert(length(dir_cell_list) <=3,'only 3 elements OR less allowed');
            
            % check dir_ang_list
            assert(length(dir_ang_list) <=3,'only 3 elements OR less allowed');
            
            % set the NUM_ROTATIONS property
            N = length(dir_ang_list);
            OBJ.num_rotations = N;
            
            % set the angle properties
            if(OBJ.ang_units=='SYM')
                tmp_ang_list = sym([0;0;0]);
            else
               tmp_ang_list       = [0;0;0];
            end
            tmp_ang_list(1:N) = dir_ang_list;
            dir_ang_list      = tmp_ang_list;
            
            OBJ.ang_1st = dir_ang_list(1);
            OBJ.ang_2nd = dir_ang_list(2);
            OBJ.ang_3rd = dir_ang_list(3);
           
            % set the direction properties
               tmp_dir_list = {'---';'---';'---'};
               tmp_dir_list(1:N) = dir_cell_list;
            dir_cell_list = tmp_dir_list;
            
           OBJ.dir_1st = dir_cell_list{1};
           OBJ.dir_2nd = dir_cell_list{2};
           OBJ.dir_3rd = dir_cell_list{3};
        end
        %------------------------------------------------------------------ 
        function OBJ = set.dir_1st(OBJ, val)
                OBJ.dir_1st = LOC_assert_valid_dir_ID(OBJ, val, '1'); 
        end
        %------------------------------------------------------------------ 
        function OBJ = set.dir_2nd(OBJ, val)
                OBJ.dir_2nd = LOC_assert_valid_dir_ID(OBJ, val, '2'); 
        end
        %------------------------------------------------------------------         
        function OBJ = set.dir_3rd(OBJ, val)
                OBJ.dir_3rd = LOC_assert_valid_dir_ID(OBJ, val, '3'); 
        end
        %------------------------------------------------------------------
        function OBJ = set.ang_units(OBJ, val)
            validateattributes(val,{'char'},{'vector'});
            val = upper(val);
            tf_test = strcmp(val, {'DEGREES','RADIANS','SYM'});
            
            assert(any(tf_test), 'ERROR:  UNknown units string');
            
            OBJ.ang_units = categorical({val});
        end
        %------------------------------------------------------------------
        function [R1] = get_active_R1(OBJ)
            
            a = OBJ.ang_1st;
            % convert to RADIANS
            if(OBJ.ang_units == 'DEGREES')
                a = a * pi/180;
            end
            
            switch OBJ.dir_1st
                case 'D1X'
                           R1 = LOC_active_RX(a);
                case 'D1Y'
                           R1 = LOC_active_RY(a);
                case 'D1Z'
                           R1 = LOC_active_RZ(a);
                case '---'
                           R1 = LOC_RI;
                otherwise
                    error('ERROR;  something UNexpected happened');
            end %switch
        end        
        %------------------------------------------------------------------
        function [R2] = get_active_R2(OBJ)
            
            a = OBJ.ang_2nd;
            % convert to RADIANS
            if(OBJ.ang_units == 'DEGREES')
                a = a * pi/180;
            end
            
            switch OBJ.dir_2nd
                case 'D2X'
                           R2 = LOC_active_RX(a);
                case 'D2Y'
                           R2 = LOC_active_RY(a);
                case 'D2Z'
                           R2 = LOC_active_RZ(a);
                case '---'
                           R2 = LOC_RI;
                otherwise
                    error('ERROR;  something UNexpected happened');
            end %switch
        end        
        %------------------------------------------------------------------
        function [R3] = get_active_R3(OBJ)
            
            a = OBJ.ang_3rd;
            % convert to RADIANS
            if(OBJ.ang_units == 'DEGREES')
                a = a * pi/180;
            end
            
            switch OBJ.dir_3rd
                case 'D3X'
                           R3 = LOC_active_RX(a);
                case 'D3Y'
                           R3 = LOC_active_RY(a);
                case 'D3Z'
                           R3 = LOC_active_RZ(a);
                case '---'
                           R3 = LOC_RI;
                otherwise
                    error('ERROR;  something UNexpected happened');
            end %switch
        end        
        %------------------------------------------------------------------
        function [aR1R2] = get_active_R1R2(OBJ)
              aR1   = get_active_R1(OBJ);
              aR2   = get_active_R2(OBJ);
            
            aR1R2   = aR1 * aR2;
        end
        %------------------------------------------------------------------
        function [aR1R2R3] = get_active_R1R2R3(OBJ)
              aR3   = get_active_R3(OBJ);
              aR2   = get_active_R2(OBJ);
              aR1   = get_active_R1(OBJ);
            
            aR1R2R3 = aR1 * aR2 * aR3;
        end
        %------------------------------------------------------------------
        function aR = get_active_R(OBJ)
              aR3   = get_active_R3(OBJ);
              aR2   = get_active_R2(OBJ);
              aR1   = get_active_R1(OBJ);
            
              aR    = aR1 * aR2 * aR3;
        end
        %------------------------------------------------------------------
        function new_v = apply_1st_rot(OBJ,v)
            new_v = apply_active_R1(OBJ, v);
        end
        %------------------------------------------------------------------        
        function new_v = apply_1st_and_2nd_rot(OBJ,v)
            new_v = apply_active_R1R2(OBJ, v);
        end 
        %------------------------------------------------------------------        
        function new_v = apply_1st_and_2nd_and_3rd_rot(OBJ,v)
            new_v = apply_active_R1R2R3(OBJ, v);
        end 
        %------------------------------------------------------------------
        function new_v = apply_active_R1(OBJ, v)
            validateattributes(v, {'double', 'sym'}, {'nrows',3})     
            aR1    = get_active_R1(OBJ);
            new_v  = aR1 * v;
        end
        %------------------------------------------------------------------
        function new_v = apply_active_R1R2(OBJ, v)
            validateattributes(v, {'double', 'sym'}, {'nrows',3})        
            aR1R2  = get_active_R1R2(OBJ);
            new_v = aR1R2 * v;
        end
        %------------------------------------------------------------------
        function new_v = apply_active_R1R2R3(OBJ, v)
            validateattributes(v, {'double', 'sym'}, {'nrows',3})       
            aR1R2R3  = get_active_R1R2R3(OBJ);
            new_v = aR1R2R3 * v;
        end
        %------------------------------------------------------------------
        function str = get_description(OBJ)
            
            tmp1_str = [char(OBJ.dir_1st),'=',num2str(OBJ.ang_1st,'%.2f')];
            tmp2_str = [char(OBJ.dir_2nd),'=',num2str(OBJ.ang_2nd,'%.2f')];
            tmp3_str = [char(OBJ.dir_3rd),'=',num2str(OBJ.ang_3rd,'%.2f')];
            tmp4_str = [char(OBJ.ang_units)];
            
            str = sprintf('%s, %s, %s, [%s]', ...
                tmp1_str, tmp2_str, tmp3_str, tmp4_str);
        end
        
    end % methods
    
end % classdef
%_#########################################################################
% END of CLASSDEF
%_#########################################################################
function out_val = LOC_assert_valid_dir_ID(OBJ, val, num_order_ID_str)
   validateattributes(val,{'char'},{'numel',3});
   val = upper(val); 
   
   % first test
   tf_test = (val=='D1X' | val=='D1Y' | val=='D1Z' | ...
              val=='D2X' | val=='D2Y' | val=='D2Z' | ...
              val=='D3X' | val=='D3Y' | val=='D3Z' | ...
              val=='---'  );
   tf_test = all(tf_test);       
   assert(tf_test,'ERROR: you can only specify D{123}{XYZ} or ---, eg: D1X or D3Z');
   
   % next test
   tf_test = val(2)==num_order_ID_str | val(2)=='-';
   assert(tf_test,'ERROR: inconsistant order ID');
   
   out_val = categorical({val});
end
%--------------------------------------------------------------------------
function aRZ = LOC_active_RZ(a)
     RZ  = LOC_RZ(a);
    aRZ  = RZ.';
end
%--------------------------------------------------------------------------
function aRY = LOC_active_RY(a)
     RY  = LOC_RY(a);
    aRY  = RY.';
end
%--------------------------------------------------------------------------
function aRX = LOC_active_RX(a)
     RX  = LOC_RX(a);
    aRX  = RX.';
end
%--------------------------------------------------------------------------
function RZ = LOC_RZ(a) % PASSIVE

    RZ = [  cos(a),  sin(a),  0;
           -sin(a),  cos(a),  0;
                 0,       0,  1;
         ];
end
%--------------------------------------------------------------------------
function RY = LOC_RY(a) % PASSIVE

    RY = [  cos(a),      0,   -sin(a);
                 0,      1,         0;
            sin(a),      0,    cos(a);
         ];
end
%--------------------------------------------------------------------------
function RX = LOC_RX(a) % PASSIVE

    RX = [       1,          0,         0;
                 0,     cos(a),    sin(a);
                 0,    -sin(a),    cos(a);
         ];
end
%--------------------------------------------------------------------------
function RI = LOC_RI()

    RI = [  1,   0,    0;
            0,   1,    0;
            0,   0,    1;
         ];
end
%--------------------------------------------------------------------------
