


for i=1:9
    %close all
    T1 = [];
    T2 = [];
    N1 = round(rand(1)*100);
    N2 = round(rand(1)*100);
    
    if N1 == 0
         N1 = round(rand(1)*100);
    end
    if N2 == 0
         N2 = round(rand(1)*100);
    end
    T1=getrow(T,T.neuron==N1);
    T2=getrow(T,T.neuron==N2);
    
    
    createElipse(i, T1, T2, numStim)
        
end