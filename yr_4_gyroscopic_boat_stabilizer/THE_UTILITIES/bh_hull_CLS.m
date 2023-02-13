classdef bh_hull_CLS
    properties
        Lix;
        Liy;
        tx;
        ty;
    end
%**************************************************************************
    methods
        function OBJ = bh_hull_CLS(Lix,Liy,tx,ty)
            arguments
                Lix (1,1) double = 1
                Liy (1,1) double = 1
                tx (1,1) double = 0.1
                ty (1,1) double = 0.1
            end
            
            OBJ.Lix = Lix;
            OBJ.Liy = Liy;
            OBJ.tx = tx;
            OBJ.ty = ty;
        end
        %------------------------------------------------------------------
        function mat = get_xy_outline_mat(obj)
           mat = [ ...
                       0,           0;
                  2*(obj.Lix + obj.tx),           0;
                  2*(obj.Lix + obj.tx),  (obj.Liy + obj.ty);
                  (2*obj.Lix + obj.tx),  (obj.Liy + obj.ty);
                  (2*obj.Lix + obj.tx),  obj.ty;
                  obj.tx, obj.ty;
                  obj.tx, (obj.ty + obj.Liy);
                  0, (obj.ty + obj.Liy);
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

