function IPM=rsa_squareIPM(IPM_vec)
% converts set of CMs (stacked along the 3rd dimension)
% to lower-triangular form (set of row vectors)
if isstruct(IPM_vec)
   N=length(IPM)
   [n,m]=size(IPM_vec(1).IPM);
    if (n~=1)
        error('IPM need to be row-vectors'); 
    end; 
    K = floor(sqrt(m*2));
    if (K*(K-1)/2+K~=m) 
        error('bad vector size'); 
    end; 
    indx=tril(true(K),0);
    IPM=IPM_vec;
    for i=1:nIPM
        A            = zeros(K);                      % make matrix 
        A(indx>0)    = IPM_vec(1,:,i);                % Assign the lower triag
        Atrans       = A';                % Now make symmetric           
        A(A==0)      = Atrans(A==0);    
        IPM(i).IPM   = A; 
    end
else
    % bare
    [n,m,nIPM]=size(IPM_vec);
    if (n~=1)
        error('IPM need to be row-vectors'); 
    end; 
    K = floor(sqrt(m*2));
    if (K*(K-1)/2+K~=m) 
        error('bad vector size'); 
    end; 
    indx=tril(true(K),0);
    IPM=[];
    for i=1:nIPM
        A            = zeros(K);                      % make matrix 
        A(indx>0)    = IPM_vec(1,:,i);                % Assign the lower triag
        Atrans       = A';                % Now make symmetric           
        A(A==0)      = Atrans(A==0);    
        IPM(:,:,i)   = A; 
    end
end