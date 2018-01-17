 function m=metric(obj,X)
            m=4./(1-X(1:2:end).^2-X(2:2:end).^2).^2;
        end