function computeAs(obj)
    if isempty(obj.V2A)
        obj.computeV2A();
    end
    fprintf('==== Computing differentials ===\n');
    tid=tic;
    %obj.As = permute(reshape(obj.V2A*obj.flat_V,2,[],2),[1 3 2]);

    obj.As = reshape(obj.V2A*obj.flat_V(:),2,2,[]);
    fprintf('done computing, %f seconds\n',toc(tid));
end