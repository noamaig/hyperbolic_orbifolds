
function  visualize( obj,varargin)
parser = inputParser;
parser.KeepUnmatched=true;
parser.addOptional('dim',2,@isnumeric);
parser.addOptional('V',[]);
parser.addOptional('landmarks',true,@islogical);
parser.addOptional('boundary',true,@islogical);
parser.addOptional('boundary_width',3,@isnumeric);
parser.addOptional('landmark_radius',1,@isnumeric);
parser.addOptional('triangle_edges',false,@islogical);
parser.addOptional('gradient',false,@islogical);
parser.addOptional('numbers',false,@islogical);
parser.addOptional('face_color','vcol');
parser.addOptional('shadow',true,@islogical);
parser.addOptional('edge_alpha',0.05,@isnumeric);
parser.addOptional('face_alpha',1);
parser.addOptional('h_model','poincare');
parser.addOptional('edge_width',0.2,@isnumeric);
parser.addOptional('circle_flips',false,@islogical);
parser.addOptional('consistent_cut_colors',false,@islogical);
parser.parse(varargin{:});

hold on;
if ~isempty(parser.Results.V)
    curY=parser.Results.V;
elseif parser.Results.dim==2
    curY=obj.flat_V;
    if strcmp(parser.Results.h_model,'klein')
        curY=poincare_to_klein(obj.flat_V);
    end
else
    curY=obj.M_cut.V;
end
if parser.Results.gradient
    X=obj.flat_V';
    X=X(:);
    [~,g]=obj.solver.objective_and_grad(X);
    g=reshape(g,2,round(length(g)/2))';
    g=g/max(abs(g(:)));
    quiver(obj.flat_V(:,1),obj.flat_V(:,2),g(:,1),g(:,2));
    return;
end
T=obj.flat_T;



    
    col=computeVCol(obj,false,varargin{:});
    %     if parser.Results.shadow
    %         col=col.*[ncol ncol ncol];
    %     end
    p=patch('vertices',curY,'faces',T,'FaceVertexCData',col,'FaceColor','interp','edgecolor','none');%,'linewidth',0.1,'edgealpha',0.05);

    set(p,'edgecolor','k','LineWidth',0.1,'edgealpha',0.2);


    [a,b]=meshgrid(1:length(obj.uncut_cone_inds), 1:length(obj.uncut_cone_inds));
    pairs=[a(:) b(:)];
    pairs=sort(pairs,2);
    pairs(pairs(:,1)==pairs(:,2),:)=[];
    pairs=unique(pairs,'rows');
    
    
    cols=linspecer(length(pairs));
    line_colors=[];
    for i=1:length(obj.M_cut.pathPairs)
        p=obj.M_cut.pathPairs{i}([1 end],1);
        p=obj.M_cut.New2Old(p);
        for j=1:2
            p(j)=find(p(j)==obj.uncut_cone_inds);
        end
        p=sort(p);
        ind=find(ismember(pairs,p,'rows'));
        line_colors(i,:)=cols(ind,:);
    end



    if obj.isdisc
        line_colors=linspecer(length(obj.reflection_paths));
        for i=1:length(obj.reflection_paths)
            p=obj.reflection_paths{i};
            %scatter(curY(p(1,1),1),curY(p(1,1),2),80,'filled');
            c=line_colors(i,:);%hsv2rgb([i/length(path_pairs),1,1]);
            
            if size(curY,2)==2
                
                line(curY(p,1),curY(p,2),'linewidth',parser.Results.boundary_width,'Color',c);
                
                
            else
                
                line(curY(p,1),curY(p,2),curY(p,3),'linewidth',3,'Color',c);
            end
        end
    else
        for i=1:length(obj.M_cut.pathPairs)
            p=obj.M_cut.pathPairs{i};
            %scatter(curY(p(1,1),1),curY(p(1,1),2),80,'filled');
            c=line_colors(i,:);%hsv2rgb([i/length(path_pairs),1,1]);
            
            if size(curY,2)==2
                
                line(curY(p(:,1),1),curY(p(:,1),2),'linewidth',parser.Results.boundary_width,'Color',c);
                line(curY(p(:,2),1),curY(p(:,2),2),'linewidth',parser.Results.boundary_width,'Color',c);
                cwidth=0.05;
                
                
            else
                
                line(curY(p(:,1),1),curY(p(:,1),2),curY(p(:,1),3),'linewidth',3,'Color',c);
                line(curY(p(:,2),1),curY(p(:,2),2),curY(p(:,1),3),'linewidth',3,'Color',c);
            end
        end
    end

% caxis([1 c_lim])
axis equal
%saveas(h,sprintf('%d.png',i),'png');
    for i=1:length(obj.uncut_cone_inds)
        if ~obj.isdisc
            curInds=obj.M_cut.Old2New{obj.uncut_cone_inds(i)};
        else
            curInds=obj.cut_cone_inds(i);
        end
        for j=1:length(curInds)
            p=curInds(j);
            cp=curY(p,:);
                   
            obj.drawLandmark(cp,obj.LM_colors(i,:),parser.Results.landmark_radius);
            if parser.Results.numbers
                text(cp(1),cp(2),cp(3),num2str(i));
            end
        end
        
        
        
    end

hold off

end
