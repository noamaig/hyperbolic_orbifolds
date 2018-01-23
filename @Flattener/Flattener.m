classdef Flattener < handle
    %Main class for embedding
    
    properties
        %the optimizer
        solver;
        %the laplacian matrix used for the weights
        L;
        %the original mesh
        M_orig;
        %indices of the cones on the original mesh
        uncut_cone_inds;
        %the indices of the vertices that were cones on the cut mesh
        cut_cone_inds;
        %positions of vertices in embedding
        flat_V;
        %the triangles of the cut mesh
        flat_T;
        %the cut mesh
        M_cut;
        %array, =1 where flipped
        flipped;
        %matrix taking vertices to their differentials
        V2A;
        %the differentials
        As;
        %determinant of differentials
        dets;
        %frobenius norm of differentials
        frobenius;
        %min sing val
        smin;
        %max sing val
        smax;
        %areas of differentials
        areas;
        
        cut_colors={[1 0 0],[0 1 0],[0 0 1],[0 1 1], [1 0 1], [1 1 0]};
        LM_colors={[1 0 0],[0 1 0],[0 0 1],[0 1 1], [1 0 1], [1 1 0]};
        orgBoundary=[];
        isDisc=false;
        lambda=1;
        side=true;
        P;
        M;
        times=[];
        LM_dist;
        mesh_name;
        reorder_cones;
        uv;
        
        reflection_paths;
        reflection_normals;
        reflections_shifts;
        reflection_normals_per_vertex=[0 0];
        reflection_shifts_per_vertex=0;
        isdisc=false;
        cone_permutation;
        auto_compute_distortion=false;
    end
    
    methods
        function obj=Flattener(V,T,uncut_cone_inds,varargin)
            % V - vertices
            % T - tris
            % uncut_cone_inds - indices of the cones
            % + other possible arguments
            parser=inputParser();
            parser.addOptional('M_cut',[]);
            parser.addOptional('isdisc',false,@islogical);
            parser.parse(varargin{:});
            if size(uncut_cone_inds,1)==1
                uncut_cone_inds=uncut_cone_inds';
            end
            assert(length(uncut_cone_inds)>=5,'Hyperbolic orbifolds have at least 5 cones!');
            %is this a disk orbifold?
            obj.isdisc=parser.Results.isdisc;
            obj.M_orig=[];
            if size(V,1)~=3 && size(V,1)~=2
                V=V';
            end
            if size(T,1)~=3
                T=T';
            end
            obj.M_orig.V=V;
            obj.M_orig.F=T;
            
            if ~obj.isdisc
                [ obj.P,obj.M ] = hyperbolic_polygon2(length(uncut_cone_inds ));
            else
                
                obj.P=hyperbolic_polygon2_disk(length(uncut_cone_inds )*2-2);
                obj.M_cut=[];
                obj.M_cut.T=obj.M_orig.F';
                obj.M_cut.V=obj.M_orig.V';
                obj.M_cut.pathPairs={};
                TR=triangulation(obj.M_orig.F',obj.M_orig.V');
                t=(TR.freeBoundary());
                t=t(end:-1:1,1);
                VV=inf(size(obj.M_orig.V));
                VV(:,t)=obj.M_orig.V(:,t);
                obj.uncut_cone_inds=knn(VV',obj.M_orig.V(:,uncut_cone_inds)');
                uncut_cone_inds=obj.uncut_cone_inds;
                assert(all(ismember(uncut_cone_inds,t)));
                obj.orgBoundary=t;
                ind=find(t==uncut_cone_inds(1));
                t=[t(ind:end); t(1:ind-1)];
                start=t(1);
                for i=2:length(uncut_cone_inds)
                    ind=find(t==uncut_cone_inds(i));
                    tt=t(1:ind);
                    obj.reflection_paths{i-1}=tt;
                    t=t(ind:end);
                end
                obj.reflection_paths{end+1}=[t; start];
                obj.cut_cone_inds=uncut_cone_inds;
                obj.flat_T=obj.M_orig.F';
                obj.M_cut.New2Old=1:length(obj.M_orig.V);
                obj.M_cut.Old2New=num2cell(1:length(obj.M_orig.V));
            end
            %             assert(length(uncut_cone_inds)==length(singularities)+1);
            
            obj.uncut_cone_inds=uncut_cone_inds;
            if ~isempty(parser.Results.M_cut)
                obj.setCutMesh(parser.Results.M_cut);
            end
            TR=triangulation(obj.M_orig.F',obj.M_orig.V');
            t=(TR.freeBoundary());
            if ~isempty(t)
                obj.orgBoundary=t(:,1);
            end
            obj.isDisc=~isempty(obj.orgBoundary);
            obj.cut_colors=linspecer(length(obj.M));
            a=linspecer(length(uncut_cone_inds));
            obj.LM_colors=a(end:-1:1,:);
            if obj.isdisc
                obj.cut_colors=linspecer(length(uncut_cone_inds));
                a=linspecer(length(uncut_cone_inds));
                obj.LM_colors=a(end:-1:1,:);
            end
        end
        
        
        function V=toHyperboloid(obj)
            %return vertex positions on hyperboloid model
            n=sum(obj.flat_V.^2,2);
            V=[1+n 2*obj.flat_V];
            for i=1:3
                V(:,i)=V(:,i)./(1-n);
            end
        end
        function V=toKlein(obj)
            %return vertex positions on Klein model
            n=1+sum(obj.flat_V.^2,2);
            V=2*[obj.flat_V(:,1)./n obj.flat_V(:,2)./n];
        end
        
        
        
        
        
        
        
        function disk_orbifold(obj)
            %embed to disk orbifold
            tid=tic;
            tr=triangulation(obj.M_cut.T,obj.M_cut.V);
            binds=tr.freeBoundary();
            clear tr;
            binds=binds([end:-1:1],1);
            VV=inf(size(obj.M_orig.V));
            VV(:,binds)=obj.M_orig.V(:,binds);
            obj.uncut_cone_inds=knn(VV',obj.M_orig.V(:,obj.uncut_cone_inds)');
            obj.cut_cone_inds=obj.uncut_cone_inds;
            assert(all(ismember(obj.uncut_cone_inds,binds)));
            
            obj.init();
            cone_inds=obj.cut_cone_inds;
            [ Porig ] = hyperbolic_polygon2_disk(  length(cone_inds));
            P=poincare_to_klein(Porig);
            
            ind=find(binds==cone_inds(1));
            binds=[binds(ind:end); binds(1:ind-1)];
            ind2=find(binds==cone_inds(2));
            ind3=find(binds==cone_inds(3));
            if ind3<ind2
                error('cones are ordered in reverse to boundary - should be CW, not CCW');
                %                cone_inds=cone_inds(end:-1:1);
                %                obj.uncut_cone_inds=cone_inds;
                %
                %                obj.cut_cone_inds=obj.uncut_cone_inds;
                %                cone_inds=obj.uncut_cone_inds;
                binds=binds(end:-1:1);
                ind=find(binds==cone_inds(1));
                binds=[binds(ind:end); binds(1:ind-1)];
                obj.M_cut.T=obj.M_cut.T(:,[1 3 2]);
                obj.flat_T=obj.M_cut.T;
                obj.M_orig.T=obj.M_cut.T';
                %                 binds=[binds;binds(1)];
            end
            binds=[binds;binds(1)];
            
            paths={};
            ns=[];
            as=[];
            for i=2:length(cone_inds)
                ind=find(binds==cone_inds(i));
                paths{end+1}=binds(1:ind);
                binds(1:ind-1)=[];
                n=P(i,:)-P(i-1,:);
                n=n/norm(n);
                n=n([2 1]);
                n(1)=-n(1);
                a=n*P(i,:)';
                ns=[ns;n];
                as=[as;a];
            end
            paths{end+1}=binds;
            n=P(1,:)-P(end,:);
            n=n/norm(n);
            n=n([2 1]);
            n(1)=-n(1);
            a=n*P(1,:)';
            ns=[ns;n];
            as=[as;a];
            obj.reflection_paths=paths;
            obj.reflection_normals=ns;
            obj.reflections_shifts=as;
            for i=1:length(paths)
                p=paths{i};
                obj.reflection_normals_per_vertex(p,:)=repmat(ns(i,:),length(p),1);
                obj.reflection_shifts_per_vertex(p)=repmat(as(i),length(p),1);
            end
            
            
            obj.solver=SolverKlein(obj.L,paths,ns,as);
            obj.klein_pseudotutte();
            x0=obj.toKlein();
            %             A=[kron(obj.L,eye(2)) obj.solver.A';
            %             obj.solver.A zeros(size(obj.solver.A,1))];
            %             b=[zeros(size(obj.L,1)*2,1);obj.solver.b];
            %             x0=A\b;
            %             x0=x0(1:length(obj.L)*2);
            %              X=x0(1:2:end);
            %             Y=x0(2:2:end);
            %             newX=[X Y];
            %              obj.flat_V=klein_to_poincare(newX);
            %              obj.computeDistortion(true);
            %              close all;
            %              obj.visualize('V',newX);
            %              pause
            x0=x0';
            x0=x0(:);
            
            x=obj.solver.solve_bfgs(x0);
            X=x(1:2:end);
            Y=x(2:2:end);
            newX=[X Y];
            %                 for i=1:length(obj.M_cut.pathPairs)
            %                     p=obj.M_cut.pathPairs{i};
            %                     newX(p(:,2),:)=obj.M{i}.map(newX(p(:,1),:));
            %                 end
            
            
            obj.flat_V=klein_to_poincare(newX);
            obj.times.flatten_orbifold=toc(tid);
            
        end
        function hyperbolic_tutte(obj)
            tid=tic;
            %convexBoundary - 1. omitted or false for free boundary
            %2. 'square' for arrangment on a square
            %3. true for disc
            %show the chosen indices
            %             if isempty(obj.cut_cone_inds)
            obj.init();
            %             end
            V=obj.hyperbolic_convex_boundary();
            inds=find(~isnan(V(:,1)));
            V=klein_to_poincare(V(inds,:));
            obj.solver=bfgs_solver(obj.L,obj.M_cut.pathPairs,obj.M,inds,V);
            %             obj.solver=Solver(obj.L,obj.M_cut.pathPairs,obj.M,inds,V);
            
            fprintf('*** flattening: ');
            obj.klein_pseudotutte();
            
            %             else
            x0=obj.flat_V';
            x0=x0(:);
            x=obj.solver.solve_bfgs(x0);
            % x=obj.solver.LM_step(x0);
            
            X=x(1:2:end);
            Y=x(2:2:end);
            newX=[X Y];
            %                 for i=1:length(obj.M_cut.pathPairs)
            %                     p=obj.M_cut.pathPairs{i};
            %                     newX(p(:,2),:)=obj.M{i}.map(newX(p(:,1),:));
            %                 end
            
            
            obj.flat_V=newX;
            n=sqrt(sum(obj.flat_V.^2,2));
            inds=n>=1;
            
            if any(inds)
                error('vertices outside of poincare disc');
            end
            obj.times.hyperbolic_tutte=toc(tid);
            
        end
        
        function flatten_orbifold(obj,useshahar)
            if nargin<2
                useshahar=true;
            end
            tid=tic;
            if isempty(obj.cut_cone_inds)
                obj.init();
            end
            obj.solver=Solver(obj.L,obj.M_cut.pathPairs,obj.M,obj.cut_cone_inds,obj.P);
            
            
            fprintf('*** flattening: ');
            obj.klein_pseudotutte();
            x0=obj.flat_V';
            x0=x0(:);
           
            if useshahar
                x=obj.solver.solve_bfgs_fast(x0);
            else
                x=obj.solver.solve_bfgs(x0);
            end
            
            X=x(1:2:end);
            Y=x(2:2:end);
            newX=[X Y];
            
            obj.flat_V=newX;
            
            n=sqrt(sum(obj.flat_V.^2,2));
            inds=n>=1;
            
            if any(inds)
                error('vertices outside of poincare disc');
            end
            obj.times.flatten_orbifold=toc(tid);
            if obj.auto_compute_distortion
                obj.computeDistortion();
            end
            
        end
        function fixFlips(obj)
            x=obj.flat_V';
            x=x(:);
            for iter=1:500
                
                X=x(1:2:end);
                Y=x(2:2:end);
                newX=[X Y];
                obj.flat_V=newX;
                obj.computeDistortion(true);
                inds=unique(obj.flat_T(obj.flipped,:));
                curInds=inds;
                for r=1:4
                    iinds=find(sum(obj.solver.adj(curInds,:)));
                    inds=[inds;iinds'];
                    curInds=iinds';
                end
                inds=unique(inds);
               
                centroid=mean(newX(inds,:));
                Mo=disc_mobius_translation(centroid,0);
                M=obj.M;
                for i=1:length(M)
                    M{i}=Mo.compose(M{i});
                end
                P=Mo.map(obj.P);
                V=Mo.map(obj.flat_V);
                fixinds=setdiff(1:length(obj.flat_V),inds);
                flip_solver=Solver(obj.L,obj.M_cut.pathPairs,M,obj.cut_cone_inds,P,'default_x',V,'fix_inds',fixinds);
                %                 flip_solver.radius=1e+5;
                flip_solver.precondition=true;
                V=V';
                V=V(:);
                x=flip_solver.gd(V);
                X=[x(1:2:end) x(2:2:end)];
                X=Mo.inverse().map(X);
                X=X';
                x=X(:);
                X=x(1:2:end);
                Y=x(2:2:end);
                newX=[X Y];
                obj.flat_V(inds,:)=newX(inds,:);
                x=obj.flat_V';
                x=x(:);
                figure(500)
                clf
                patch('Faces',obj.flat_T(obj.flipped,:),'Vertices',obj.flat_V,'FaceColor','none');
                pause(0.2)
            end
        end
        function fixFlipsNew(obj)
            %this function fixes any flipped tris (these shouldn't happen
            %mathematically but may happen due to numerical errors). 
            
            
            tr=triangulation(obj.flat_T,obj.flat_V);
            binds=tr.freeBoundary();
            binds=unique(binds(:));
            for iter=1:10000
                
                x=obj.flat_V';
                x=x(:);
                
                
                X=x(1:2:end);
                Y=x(2:2:end);
                newX=[X Y];
                obj.flat_V=newX;
                
                %recalculate distortion
                obj.computeDistortion(true);
                %all verticers in flipped tris
                inds=unique(obj.flat_T(obj.flipped,:)');
                %no vertices means no flips :)
                if isempty(inds)
                    return;
                end
                
                %take the 5-ring neighbourhood of flipped vertices
                curInds=inds;
                for r=1:5
                    iinds=find(sum(obj.solver.adj(curInds,:)));
                    inds=[inds;iinds'];
                    curInds=iinds';
                end
                %remove duplicates
                inds=unique(inds);
                %don't take boundary vertices
                inds=setdiff(inds,binds);
                
                %iterate over all flipped vertices
                for i=1:length(inds)
                    ind=inds(i);
                    %mobius trans that takes cur vertex to origin
                    Mo=disc_mobius_translation(newX(ind,:),0);
                    %map all vertices w.r.t. to the mobius
                    tempX=Mo.map(newX);
                    w=obj.solver.Wmat(ind,:);
                    w(ind)=0;
                    %                         w=double(w~=0);
                    w=w/sum(w);
                    xx=1*w*tempX;
                    newX(ind,:)=Mo.inverse().map(xx);
                    
                    obj.flat_V=newX;
                    
                end
                %                     end
                %                 figure(500);
                %                 clf
                %                 %             patch('Faces',obj.flat_T(obj.flipped,:),'Vertices',obj.flat_V,'FaceColor','none');
                %                 obj.visualize('landmarks',false);
                %                 axis(ax);
                %                 pause(0.2);
            end
        end
        
        function fixFlipsDisk(obj)
            tr=triangulation(obj.flat_T,obj.flat_V);
            binds=tr.freeBoundary();
            binds=unique(binds(:));
            for iter=1:10000
                
                x=obj.flat_V';
                x=x(:);
                
                
                X=x(1:2:end);
                Y=x(2:2:end);
                newX=[X Y];
                obj.flat_V=newX;
                obj.computeDistortion(true);
                inds=unique(obj.flat_T(obj.flipped,:)');
                curInds=inds;
                if isempty(inds)
                    return;
                end
                for r=1:1
                    iinds=find(sum(obj.solver.adj(curInds,:)));
                    inds=[inds;iinds'];
                    curInds=iinds';
                end
                
                inds=unique(inds);
                inds=setdiff(inds,obj.uncut_cone_inds);
                
                %                     for innerinter=1:20
                for i=1:length(inds)
                    
                    ind=inds(i);
                    
                    Mo=disc_mobius_translation(newX(ind,:),0);
                    tempX=Mo.map(newX);
                    w=obj.solver.Wmat(ind,:);
                    w(ind)=0;
                    %                         w=double(w~=0);
                    w=w/sum(w);
                    xx=1*w*tempX;
                    isboundary=false;
                    if ismember(ind,binds)
                        
                        isboundary=true;
                        
                        found=0;
                        for j=1:length(obj.reflection_paths)
                            if ismember(ind,obj.reflection_paths{j})
                                found=j;
                                break;
                            end
                        end
                        assert(found~=0);
                        s=tempX(obj.reflection_paths{found}(1),:);
                        e=tempX(obj.reflection_paths{found}(end),:);
                        n=normr(e-s);
                        n=n([2 1]);
                        n(1)=-n(1);
                        
                        %                         shifts=tempX*n';
                        
                        %                         tempXR=tempX-2*kron(shifts,n);
                        xxR=xx-2*(xx*n')*n;
                        
                        xx=(xx+xxR)/2;
                    end
                    newX(ind,:)=Mo.inverse().map(xx);
                    %                     if isboundary
                    %                         XXX=poincare_to_klein(newX(ind,:));
                    %                         aa=obj.reflection_shifts_per_vertex(ind)-XXX*obj.reflection_normals_per_vertex(ind,:)';
                    %                         XXX=XXX+aa*obj.reflection_normals_per_vertex(ind,:);
                    %                         newX(ind,:)=klein_to_poincare(XXX);
                    %                     end
                    obj.flat_V=newX;
                    
                    %                     figure(500);
                    %                     ax=gca;
                    %                                clf
                    %                                obj.visualize('landmarks',false);
                    %                                axis(ax);
                    %                                pause;
                    %             patch('Faces',obj.flat_T(obj.flipped,:),'Vertices',obj.flat_V,'FaceColor','none');
                    %                 obj.visualize('landmarks',false);
                    %                 axis(ax);
                    %                 pause(0.2);
                end
                %                     end
                %
            end
        end
        
        
        
        
        function V=hyperbolic_convex_boundary(obj)
            if obj.isdisc
                obj.P= hyperbolic_polygon2_disk(  length(obj.cut_cone_inds));
            end
            klein_P=poincare_to_klein(obj.P);
            V=nan(length(obj.M_cut.V),2);
            
            V(obj.cut_cone_inds,:)=klein_P;
            if ~obj.isdisc
                for j=1:length(obj.M_cut.pathPairs)
                    %take the two paths
                    path1=obj.M_cut.pathPairs{j}(:,1);
                    path2=obj.M_cut.pathPairs{j}(:,2);
                    
                    
                    %indices of previous and next point
                    
                    %if this is the first iteration, set the position of the
                    %boundary to a triangle
                    
                    x1=linspace(V(path1(1),1),V(path1(end),1),length(path1));
                    y1=linspace(V(path1(1),2),V(path1(end),2),length(path1));
                    v1=[x1' y1'];
                    V(path1,:)=v1;
                    %
                    
                    v2=klein_to_poincare(v1);
                    
                    v2=obj.M{j}.map(v2);
                    v2=poincare_to_klein(v2);
                    
                    
                    V(path2,:)=v2;
                    
                end
            else
                for j=1:length(obj.reflection_paths)
                    %take the two paths
                    path1=obj.reflection_paths{j};
                    
                    
                    
                    %indices of previous and next point
                    
                    %if this is the first iteration, set the position of the
                    %boundary to a triangle
                    
                    x1=linspace(V(path1(1),1),V(path1(end),1),length(path1));
                    y1=linspace(V(path1(1),2),V(path1(end),2),length(path1));
                    v1=[x1' y1'];
                    V(path1,:)=v1;
                end
            end
        end
        function A=landmark_distances(obj)
            if ~isempty(obj.LM_dist)
                A=obj.LM_dist;
                return;
            end
            A=zeros(length(obj.uncut_cone_inds));
            G=adjacency_edge_cost_matrix(obj.M_orig.V',obj.M_orig.F');
            if exist('graph')%check if matlab has the graph library
                G=graph(G);
                A=G.distances(obj.uncut_cone_inds,obj.uncut_cone_inds);
                return;
            end
            for i=1:length(obj.uncut_cone_inds)
                i
                for j=i+1:length(obj.uncut_cone_inds)
                    
                    A(i,j) = graphshortestpath(G,obj.uncut_cone_inds(i),obj.uncut_cone_inds(j));
                    A(j,i)=A(i,j);
                end
            end
            obj.LM_dist=A;
        end
        function sortByDistance(obj)
            A=obj.landmark_distances();
            A=A+sparse(1:length(A),1:length(A),inf);
            
            newV=1;
            v=2:length(obj.uncut_cone_inds);
            while(~isempty(v))
                [~,ind]=min(A(newV(end),v));
                newV=[newV v(ind)];
                v(ind)=[];
                
            end
            obj.uncut_cone_inds=obj.uncut_cone_inds(newV);
        end
        function w=colorWeights(obj,flat,r)
            
            D=zeros(length(obj.M_orig.V),length(obj.uncut_cone_inds));
            if flat
                %                 G=adjacency_edge_cost_matrix(obj.flat_V,obj.flat_T).^2;
                %
                E = edges(obj.flat_T);
                %
                %
                %   % compute edge norms
                edge_norms = sum(1000*(obj.flat_V(E(:,1),:)-obj.flat_V(E(:,2),:)).^2,2);
                E=obj.M_cut.New2Old(E);
                E=sort(E,2);
                [E,inds]=unique(E,'rows');
                edge_norms=edge_norms(inds);
                n = length(obj.M_orig.V);
                % build sparse adjacency matrix with non-zero entries indicated edge costs
                G = sparse([E(:,1);E(:,2)],[E(:,2);E(:,1)],repmat(edge_norms,[2 1]),n,n);
                %
                %
                for i=1:length(obj.uncut_cone_inds)
                    d=sqrt(graphshortestpath(G,obj.uncut_cone_inds(i)));
                    d=exp(-r*d.^2);
                    D(:,i)=d;
                end
            else
                G=adjacency_edge_cost_matrix(obj.M_orig.V',obj.M_orig.F').^2;
                %10 was good
                if nargin<3
                    r=25;
                end
                for i=1:length(obj.uncut_cone_inds)
                    d=sqrt(graphshortestpath(G,obj.uncut_cone_inds(i)));
                    d=exp(-r*d.^2);
                    D(:,i)=d;
                end
                %                 D=bsxfun(@times,D,1./sum(D,2).^0.8);
                w=D(obj.M_cut.New2Old,:);
            end
            
        end
        function findBestOrder(obj)
            addpath('C:\Program Files\Mosek\7\toolbox\r2013a');
            addpath(genpath('YALMIP-develop'));
            addpath(genpath('gptoolbox-master'));
            tic
            
            D1 = obj.landmark_distances();
            n = length(D1);
            [~,~,~,D2] =hyperbolic_polygon2(length(obj.uncut_cone_inds));
            
            
            
            % solve graph matching
            P = binvar(n,n,'full'); % as an integer program
            %             P = sdpvar(n,n,'full'); % relaxation
            h = norm(D1*P-P*D2,'fro');
            F = (P(1,1)==1) + (ones(1,n)*P==ones(1,n)) + (P*ones(n,1)==ones(n,1)) + (P>=0);
            a=sdpsettings();
            a.mosek.MSK_DPAR_MIO_MAX_TIME=600;
            solvesdp(F,h,a);
            obj.cone_permutation=double(P)';
            obj.uncut_cone_inds=obj.cone_permutation*obj.uncut_cone_inds;
        end
        function orderTS(obj)
            D1 = obj.landmark_distances();
            nStops=length(obj.uncut_cone_inds);
            idxs = nchoosek(1:nStops,2);
            
            
            dist = D1(sub2ind(size(D1),idxs(:,1),idxs(:,2)));
            lendist = length(dist);
            
            Aeq = spones(1:length(idxs)); % Adds up the number of trips
            beq = nStops;
            
            Aeq = [Aeq;spalloc(nStops,length(idxs),nStops*(nStops-1))]; % allocate a sparse matrix
            for ii = 1:nStops
                whichIdxs = (idxs == ii); % find the trips that include stop ii
                whichIdxs = sparse(sum(whichIdxs,2)); % include trips where ii is at either end
                Aeq(ii+1,:) = whichIdxs'; % include in the constraint matrix
            end
            beq = [beq; 2*ones(nStops,1)];
            
            intcon = 1:lendist;
            lb = zeros(lendist,1);
            ub = ones(lendist,1);
            
            opts = optimoptions('intlinprog','Display','off');
            [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);
            
            tours = detectSubtours(x_tsp,idxs);
            numtours = length(tours); % number of subtours
            fprintf('# of subtours: %d\n',numtours);
            A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
            b = [];
            while numtours > 1 % repeat until there is just one subtour
                % Add the subtour constraints
                b = [b;zeros(numtours,1)]; % allocate b
                A = [A;spalloc(numtours,lendist,nStops)]; % a guess at how many nonzeros to allocate
                for ii = 1:numtours
                    rowIdx = size(A,1)+1; % Counter for indexing
                    subTourIdx = tours{ii}; % Extract the current subtour
                    %         The next lines find all of the variables associated with the
                    %         particular subtour, then add an inequality constraint to prohibit
                    %         that subtour and all subtours that use those stops.
                    variations = nchoosek(1:length(subTourIdx),2);
                    for jj = 1:length(variations)
                        whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                            (sum(idxs==subTourIdx(variations(jj,2)),2));
                        A(rowIdx,whichVar) = 1;
                    end
                    b(rowIdx) = length(subTourIdx)-1; % One less trip than subtour stops
                end
                
                % Try to optimize again
                [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub,opts);
                
                
                % How many subtours this time?
                tours = detectSubtours(x_tsp,idxs);
                numtours = length(tours); % number of subtours
                fprintf('# of subtours: %d\n',numtours);
            end
            dist;
            tour=tours{1}';
            
            
            adjDist=sparse([idxs(:,1);idxs(:,2)],[idxs(:,2);idxs(:,1)],[dist;dist]);
            [~,m]=max(adjDist(sub2ind(size(adjDist),tour,[tour(2:end);tour(1)])));
            tour=[tour(m+1:end);tour(1:m)];
            obj.uncut_cone_inds=obj.uncut_cone_inds(tour);
            obj.reorder_cones=tour;
        end
        
    end
    
end
