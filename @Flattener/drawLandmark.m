function drawLandmark(obj,cp,col,radius)
            if size(cp,2)==2
                cwidth=0.035*radius;
                rectangle('position',[cp-cwidth/2 cwidth cwidth],'curvature',1,'facecolor',col,'edgecolor','k');
            else
                
                    scatter3(cp(1),cp(2),cp(3),300,col,'filled');
                [x,y,z] = sphere;
                x=radius*x/40;
                y=radius*y/40;
                z=radius*z/40;
                
                
                
                %                 surf(x+cp(1),y+cp(2),z+cp(3),'edgecolor','none','facecolor',col) % centered at (3,-2,0)
                
            end
        end