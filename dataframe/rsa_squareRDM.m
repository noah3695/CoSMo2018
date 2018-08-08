function RDM=rsa_squareRDM(RDM_vec)
% converts set of CMs (stacked along the 3rd dimension)
% to lower-triangular form (set of row vectors)
if isstruct(RDM_vec)
   N=length(RDM)
   [n,m]=size(RDM_vec(1).RDM);
    if (n~=1)
        error('RDM need to be row-vectors'); 
    end; 
    K = ceil(sqrt(m*2));
    if (K*(K-1)/2~=m) 
        error('bad vector size'); 
    end; 
    indx=tril(true(K),-1);
    RDM=RDM_vec;
    for i=1:nRDM
        A            = zeros(K);                      % make matrix 
        A(indx>0)    = RDM_vec(1,:,i);                % Assign the lower triag
        Atrans       = A';                % Now make symmetric           
        A(A==0)      = Atrans(A==0);    
        RDM(i).RDM   = A; 
    end
else
    % bare
    [n,m,nRDM]=size(RDM_vec);
    if (n~=1)
        error('RDM need to be row-vectors'); 
    end; 
    K = ceil(sqrt(m*2));
    if (K*(K-1)/2~=m) 
        error('bad vector size'); 
    end; 
    indx=tril(true(K),-1);
    RDM=[];
    for i=1:nRDM
        A            = zeros(K);                      % make matrix 
        A(indx>0)    = RDM_vec(1,:,i);                % Assign the lower triag
        Atrans       = A';                % Now make symmetric           
        A(A==0)      = Atrans(A==0);    
        RDM(:,:,i)   = A; 
    end
end