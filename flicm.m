function [Uout,iter] = flicm( image, U, m, cNum, maxIter, thrE )
    Uold = U;
    c = ClaCenter(image,U,cNum,m);
	sweeps = 0;
	dMax = 10.0;
    tempSum = zeros(size(image,1),size(image,2),cNum);    
    dist = sqrt([2.0; 1.0; 2.0; 1.0; 1.0; 2.0; 1.0; 2.0]);     % sapce distance for neigbour
    S = uint32(1:size(image,1)*size(image,2));
    SS = reshape(S,size(image,1),size(image,2));
    padNum = 1;
    SSpad = wextend('2D','zpd',SS,padNum);
    SSpad = SSpad(:);
    M = size(image,1) + padNum*2;
    
    while( dMax>thrE && sweeps<=maxIter )
        for k = 1:cNum
            temp = shiftdim( repmat(c(k,:),[size(image,2),1,size(image,1)]) , 2 );            
            tempSum(:,:,k) = sum( (image - temp).^2, 3 );
        end
        tempSum1 = reshape(tempSum,size(tempSum,1)*size(tempSum,2),size(tempSum,3));
        sSum = zeros(1,cNum);
        Uold = reshape(Uold,size(U,1)*size(U,2),size(U,3));
        U = zeros(numel(S),cNum);
        for i = 1:numel(S)
%             claculate degree of membership for every pixel
            idx = find_8Neighbur(SSpad,i,M);
            for k = 1:cNum
                if(nnz(idx<=0))
                    idx1 = idx(idx>0);
                    sSum(k) = sum( ...
                    1.0./(1+dist(idx>0)).*(( 1 - Uold(idx1,k) ).^m).*tempSum1(idx1,k)...
                    ) + tempSum1(i,k);
                else
                    sSum(k) = sum( ...
                        1.0./(1+dist).*(( 1 - Uold(idx,k) ).^m).*tempSum1(idx,k)...
                        ) + tempSum1(i,k);
                end
            end
            meshSum = meshgrid(sSum,sSum);
            transMeshSum = meshSum';
            U(i,:) = 1.0./( sum( (meshSum./transMeshSum).^(1.0/(m-1.0)) ) );
        end
        U = reshape(U,size(image,1),size(image,2),cNum);
        c = ClaCenter(image,U,cNum,m);
        
        dMax = max(abs( Uold(:)-U(:) ));
        Uold = U;
        Uout = U;
        sweeps = sweeps + 1 ;
        disp(['Iteration ', num2str(sweeps), ': Gap ', num2str(dMax)]);
    end
    iter = sweeps;
end

function[c] = ClaCenter(image,U,cNum,m)
% claculate class center
    c = zeros(cNum,size(image,3));
    sSum = sum(sum(U.^m));
    for i = 1:cNum
        U_m = repmat(U(:,:,i).^m,[1,1,size(image,3)]);
        c(i,:) = sum(sum(U_m.*image))./sSum(i);
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
