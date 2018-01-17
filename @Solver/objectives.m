 function o=objectives(obj,X)
            if size(X,2)==1
                X=[X(1:2:end) X(2:2:end)];
            end
            
            interior=obj.karcher_objective(obj.Iinterior,obj.Jinterior,X);
            right=obj.karcher_objective(obj.Iright,obj.Jright,X);
            right_to_right=obj.karcher_objective(obj.IRightToRight,obj.JRightToRight,X);
            o=[interior;right;right_to_right];
        end