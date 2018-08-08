function IPM=rsa_squareIPMfull(IPM_vec)
% converts set of CMs (stacked along the 3rd dimension)
% to lower-triangular form (set of row vectors)
if isstruct(IPM_vec)
   N=length(IPM)
   [n,m]=size(IPM_vec(1).IPM);
    if (n~=1)
        error('IPM need to be row-vectors'); 
    end; 
    K = floor(sqrt(m));
    if (abs(K-sqrt(m))>eps) 
        error('bad vector size'); 
    end; 
    IPM=IPM_vec;
    for i=1:nIPM
        IPM(i).IPM   = reshape(IPM(i).IPM,K,K); 
    end
else
    % bare
    [n,m,nIPM]=size(IPM_vec);
    if (n~=1)
        error('IPM need to be row-vectors'); 
    end; 
    K = floor(sqrt(m));
    if (abs(K-sqrt(m))>eps) 
        error('bad vector size'); 
    end; 
    IPM=[];
    for i=1:nIPM
        IPM(:,:,i)   = reshape(IPM_vec(:,:,i),K,K); 
    end
end