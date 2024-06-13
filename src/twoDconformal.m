function [ X, RMSEcnt, xycnt, xyunk ] = twoDconformal( cntfidu,cntarbi,unkarbi )
%
%

if nargin == 2 || nargin ==3
    
    Numberofcnt = size(cntfidu,1);
    Acnt = zeros(2*Numberofcnt,4);
    Lcnt = zeros(2*Numberofcnt,1);
    
    i = 1:2:2*Numberofcnt;
    j = 2:2:2*Numberofcnt;
    
    Acnt(i,1) = +cntarbi(:,1);Acnt(i,2) = -cntarbi(:,2);Acnt(i,3) = 1;
    Acnt(j,1) = +cntarbi(:,2);Acnt(j,2) = +cntarbi(:,1);Acnt(j,4) = 1;
    
    Lcnt(i,1) = cntfidu(:,1);Lcnt(j,1) = cntfidu(:,2);
    
    X = (Acnt'*Acnt)\(Acnt'*Lcnt);
    LLcnt = Acnt*X;
    xycnt = reshape(LLcnt',2,Numberofcnt)';
    RMSEcnt = sqrt(sum((Lcnt - LLcnt).^2)/(Numberofcnt-1));
    
    if nargin == 3
        
        Numberofunk = size(unkarbi,1);
        Aunk = zeros(2*Numberofunk,4);
        
        i = 1:2:2*Numberofunk;
        j = 2:2:2*Numberofunk;
        
        Aunk(i,1) = +unkarbi(:,1);Aunk(i,2) = -unkarbi(:,2);Aunk(i,3) = 1;
        Aunk(j,1) = +unkarbi(:,2);Aunk(j,2) = +unkarbi(:,1);Aunk(j,4) = 1;
        
        LLunk = Aunk*X;
        xyunk = reshape(LLunk',2,Numberofunk)';
        
    end
end

end



