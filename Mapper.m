classdef Mapper < handle
    % Responsible for computing a map from two flattenings
    
    properties
        %array of parameterizations
        flatteners;
        %flag of whether to (==true) run TSP on first mesh and then use its landmark
        %order for the cutting, or (==false) use given order.
        computeBestOrder=true;
        %the computed map
        map;
    end
    
    methods
        function obj= Mapper(flatteners,useGivenOrder)
            % input: array of flattener objects (the parameterizations). If
            % they haven't been computed yet the mapper will flatten them
            % itself.
            % Optinal flag (default=true) as to whether use a
            % travelling-salesman approximation in order to compute the
            % best cutting order of the vertices (this is a heuristic, see
            % paper)
            obj.flatteners=flatteners;
            if nargin>1
                obj.computeBestOrder=useGivenOrder;
            end
        end
        function printout(obj,msg)
            %helper function for printing messages
            fprintf('*** Mapper: %s\n',msg);
        end
        function computeMap(obj)
            %compte the map between two objects
            obj.flatten();
            for i=2:length(obj.flatteners)
                obj.lift(1,i);
            end
        end
        function flatten(obj)
            %embed all the meshes
            for i=1:length(obj.flatteners)
                obj.printout(sprintf('flattening %d/%d',i,length(obj.flatteners)));
                if obj.computeBestOrder
                    if i==1
                        obj.flatteners{i}.orderTS();
                    else
                        obj.flatteners{i}.uncut_cone_inds=obj.flatteners{i}.uncut_cone_inds(obj.flatteners{1}.reorder_cones);
                    end
                end
                obj.flatteners{i}.flatten_orbifold();
                obj.flatteners{i}.fixFlipsNew();
            end
        end
        function lift(obj,s,t)
            %compute the surface map between all objects
            obj.map=UncutSurfMap(obj.flatteners);
            obj.map.compute(s,t);
            
        end
    end
    
end

