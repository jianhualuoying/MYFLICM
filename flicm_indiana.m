function [Uout,iter] = flicm_indiana( image, label, U, m, cNum, maxIter, thrE )

    idxNon0 = find(label ~= 0);
    img = image(idxNon0,:); 
    Uold = U;
    c = ClaCenter(img,U,cNum,m);
	sweeps = 0;
	dMax = 10.0;
    DIST = sqrt([2.0; 1.0; 2.0; 1.0; 1.0; 2.0; 1.0; 2.0]);     % sapce distance for neigbour
    
    S = uint32(1:size(image,1));
    SS = reshape(S,145,145);
    padNum = 1;
    SSpad = wextend('2D','zpd',SS,padNum);
    clear SS S;
    SSpad = SSpad(:);
    M = 145 + padNum*2;
    
    while( dMax>thrE && sweeps<=maxIter )
        
        temp = permute(c,[3,2,1]);
        tempSum = squeeze( sum( (bsxfun(@minus,image,temp)).^2, 2 ) );
        clear temp;
        
%       one interpretaion for the sentences above        
%         for k = 1:cNum
%             temp = repmat(c(k,:),[size(image,1),1]);            
%             tempSum(:,k) = sum( (image - temp).^2, 2 );
%         end

%       another interpretaion for the sentences above    
%         temp = permute(c,[3,2,1]);                              % 1D - none  2D - feature dimension 3D - class dimension
%         temp = repmat(temp,[1,size(image,1)*size(image,2),1]);  % 1D - sample dimension  2D - feature dimension 3D - class dimension
%         tempImage = repmat(image,[1,1,size(c,1)]);              % 1D - sample dimension  2D - feature dimension 3D - class dimension
%         tempSum = squeeze( sum( (tempImage - temp).^2, 2 ) );   % 1D - sample dimension  2D - class dimension
        
        sSum = zeros(1,cNum);
        row = size(img,1);
        U = zeros(row,cNum);
        for i = 1:row
            idx_i = idxNon0(i);
            idxNeighbor = find_8Neighbur(SSpad,idx_i,M);
            dist = DIST;
            if(nnz(idxNeighbor<=0))
                dist = DIST(idxNeighbor>0);
                idxNeighbor = idxNeighbor(idxNeighbor>0);
            end
            binary = ismember(idxNeighbor,idxNon0);
            dist = dist(binary);
            idxNeighbor = idxNeighbor(binary);
            idxNeighbor = find(ismember(idxNon0,idxNeighbor));
            for k = 1:cNum
                    sSum(k) = sum( ...
                    1.0./(1+dist).*(( 1 - Uold(idxNeighbor,k) ).^m).*tempSum(idxNeighbor,k)...
                    ) + tempSum(i,k);
            end
            meshSum = meshgrid(sSum,sSum);
            transMeshSum = meshSum';
            U(i,:) = 1.0./( sum( (meshSum./transMeshSum).^(1.0/(m-1.0)) ) );
        end
        
        c = ClaCenter(img,U,cNum,m);
        
        dMax = max(abs( Uold(:)-U(:) ));
        Uold = U;
        Uout = U;
        sweeps = sweeps + 1 ;
        disp(['Iteration ', num2str(sweeps), ': Gap ', num2str(dMax)]);
    end
    iter = sweeps;
end

function[c] = ClaCenter(image,U,cNum,m)
    c = zeros(cNum,size(image,2));
    sSum = sum(U.^m);
    for i = 1:cNum
        U_m = repmat(U(:,i).^m,[1,size(image,2)]);
        c(i,:) = sum(U_m.*image)./sSum(i);
    end
end

function idxVector = find_8Neighbur(X,num,cloumnNum)
% find_8Neighbur找到9宫格近邻的索引值
% cloumnNum-列的数量
% idxVector-8个近邻的索引矢量
    coreIdx = find(X==num);
    idxVector = [coreIdx-cloumnNum-1,coreIdx-cloumnNum,coreIdx-cloumnNum+1,...
        coreIdx-1,coreIdx+1,coreIdx+cloumnNum-1,coreIdx+cloumnNum,coreIdx+cloumnNum+1];
    idxVector = X(idxVector);
end
