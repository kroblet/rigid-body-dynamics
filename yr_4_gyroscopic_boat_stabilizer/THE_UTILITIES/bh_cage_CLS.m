classdef bh_cage_CLS
    properties
        Lx;
        Ly;
        tx;
        ty;
    end
%**************************************************************************
    methods
        function OBJ = bh_cage_CLS(Lx,Ly,tx,ty)
            arguments
                Lx (1,1) double = 1
                Ly (1,1) double = 1
                tx (1,1) double = 0.1
                ty (1,1) double = 0.1
            end
            
            OBJ.Lx = Lx;
            OBJ.Ly = Ly;
            OBJ.tx = tx;
            OBJ.ty = ty;
        end
        %------------------------------------------------------------------
        function mat = get_xy_outline_mat(obj)
           
           DY = 1e-7; %a small offset for the cut 
            
           mat = [ ...
                       0,            0;
                  obj.Lx,            0;
                  obj.Lx,            obj.Ly/2;  % POINT 3
                  (obj.Lx-obj.tx),   obj.Ly/2;  % POINT 4 
                  (obj.Lx-obj.tx),   obj.ty;
                  obj.tx,            obj.ty;
                  obj.tx,            (obj.Ly-obj.ty);
                  (obj.Lx-obj.tx),   (obj.Ly-obj.ty);
                  (obj.Lx-obj.tx), DY+obj.Ly/2;  % cut point for P4
                  obj.Lx,          DY+obj.Ly/2;  % cut point for P3
                  %--------------------------
                  obj.Lx,   obj.Ly;
                  0,        obj.Ly;
                  %0,        0;
                  ];
        end
        %------------------------------------------------------------------
        function plot(obj)
            mat = get_xy_outline_mat(obj);
            figure;
            %plot(mat(:,1), mat(:,2), '-ro');
            
            % why not use the POLYSHAPE() function
            P = polyshape(mat);
            P.plot();
            
            grid('on')
            hold on
            plot(mat(:,1), mat(:,2), 'k.');
            
            [Xc, Yc] = P.centroid();
            plot(Xc,Yc,'ro', "MarkerFaceColor", 'red',"MarkerEdgeColor", "k");
            
            axis('equal');
        end
        %------------------------------------------------------------------
    end
end

