classdef bh_rot_passive_G2B_CLS
% BH_ROT_PASSIVE_B2G_CLS define PASSIVE rotation matrix for converting G_vec to B_vec
%==========================================================================
% Assume that a body undergoes 3 successive rotations in the following order:
%
%   a.) R1 occurs 1st about some local body axis
%   b.) R2 occurs 2nd about some local body axis
%   c.) R3 occurs 3rd about some local body axis
%
% We can then convert a vector defined in the G axis to it's corresponding 
% description in the B axis, using a PASSIVE rotation matrix, ie:
%
%    vB = R3 * R2 * R1 * vG
%    vB = bRg * vG
%==========================================================================
% Example USAGE:
%
% % A typical usage scenario
%   OBJ_A = bh_rot_passive_G2B_CLS({'D1Z', 'D2Y', 'D3X'}, [30 20 10],'DEGREES')
% 
%   R1 = OBJ_A.get_R1
%   R2 = OBJ_A.get_R2
%   R3 = OBJ_A.get_R3
% 
%   R2R1   = OBJ_A.get_R2R1
% 
%   R3R2R1 = OBJ_A.get_R3R2R1
%--------------------------------------------------------------------------
% % convert a vector in G, into its components in B
%   vG    = [1,0,0]';
%   vB    = OBJ_A.apply_R3R2R1(vG)
%--------------------------------------------------------------------------
% %% You can also supply SYMBOLIC angles
%   OBJ_B = bh_rot_passive_G2B_CLS({'D1Z', 'D2Y', 'D3X'}, [sym('A'), sym('B'), sym('C')], 'SYM')
% 
%   R1 = OBJ_A.get_R1
% 
%   R3R2R1 = OBJ_A.get_R3R2R1
% 
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
        function OBJ = bh_rot_passive_G2B_CLS(dir_cell_list, dir_ang_list, units_str)
            
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
            
           OBJ.dir_1st = dir_cell_list{1}; % our SET methods that test it
           OBJ.dir_2nd = dir_cell_list{2}; % our SET methods that test it
           OBJ.dir_3rd = dir_cell_list{3}; % our SET methods that test it
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
        function [R1] = get_R1(OBJ)
            
            a = OBJ.ang_1st;
            % convert to RADIANS
            if(OBJ.ang_units == 'DEGREES')
                a = a * pi/180;
            end
            
            switch OBJ.dir_1st
                case 'D1X'
                           R1 = LOC_RX(a);
                case 'D1Y'
                           R1 = LOC_RY(a);
                case 'D1Z'
                           R1 = LOC_RZ(a);
                case '---'
                           R1 = LOC_RI;
                otherwise
                    error('ERROR;  something UNexpected happened');
            end %switch
        end        
        %------------------------------------------------------------------
        function [R2] = get_R2(OBJ)
            
            a = OBJ.ang_2nd;
            % convert to RADIANS
            if(OBJ.ang_units == 'DEGREES')
                a = a * pi/180;
            end
            
            switch OBJ.dir_2nd
                case 'D2X'
                           R2 = LOC_RX(a);
                case 'D2Y'
                           R2 = LOC_RY(a);
                case 'D2Z'
                           R2 = LOC_RZ(a);
                case '---'
                           R2 = LOC_RI;
                otherwise
                    error('ERROR;  something UNexpected happened');
            end %switch
        end        
        %------------------------------------------------------------------
        function [R3] = get_R3(OBJ)
            
            a = OBJ.ang_3rd;
            % convert to RADIANS
            if(OBJ.ang_units == 'DEGREES')
                a = a * pi/180;
            end
            
            switch OBJ.dir_3rd
                case 'D3X'
                           R3 = LOC_RX(a);
                case 'D3Y'
                           R3 = LOC_RY(a);
                case 'D3Z'
                           R3 = LOC_RZ(a);
                case '---'
                           R3 = LOC_RI;
                otherwise
                    error('ERROR;  something UNexpected happened');
            end %switch
        end        
        %------------------------------------------------------------------
        function [R2R1] = get_R2R1(OBJ)
              R1   = get_R1(OBJ);
              R2   = get_R2(OBJ);
            
            R2R1   = R2 * R1;
        end
        %------------------------------------------------------------------
        function [R3R2R1] = get_R3R2R1(OBJ)
              R1   = get_R1(OBJ);
              R2   = get_R2(OBJ);
              R3   = get_R3(OBJ);
            
            R3R2R1 = R3 * R2 * R1;
        end
        %------------------------------------------------------------------
        function R = get_R(OBJ)
              R1   = get_R1(OBJ);
              R2   = get_R2(OBJ);
              R3   = get_R3(OBJ);
            
            R = R3 * R2 * R1;
        end
        %------------------------------------------------------------------
        function new_v = apply_R1(OBJ, v)
            assert(iscolumn(v), 'ERROR: give me a COLUMN');
            
            R1    = get_R1(OBJ);
            new_v = R1 * v;
        end
        %------------------------------------------------------------------
        function new_v = apply_R2R1(OBJ, v)
            assert(iscolumn(v), 'ERROR: give me a COLUMN');
            
            R2R1  = get_R2R1(OBJ);
            new_v = R2R1 * v;
        end
        %------------------------------------------------------------------
        function new_v = apply_R3R2R1(OBJ, v)
            assert(iscolumn(v), 'ERROR: give me a COLUMN');
            
            R3R2R1  = get_R3R2R1(OBJ);
            new_v = R3R2R1 * v;
        end
        %------------------------------------------------------------------
        
    end
    
end
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
function RZ = LOC_RZ(a)

    RZ = [  cos(a),  sin(a),  0;
           -sin(a),  cos(a),  0;
                 0,       0,  1;
         ];
end
%--------------------------------------------------------------------------
function RY = LOC_RY(a)

    RY = [  cos(a),      0,   -sin(a);
                 0,      1,         0;
            sin(a),      0,    cos(a);
         ];
end
%--------------------------------------------------------------------------
function RX = LOC_RX(a)

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
