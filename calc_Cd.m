function Cd = calc_Cd(data,dim,distancetype,r)

    L=length(data(:,1))-dim+1;
    X=[];
    for i=1:dim
        X=[X data(i:i+L-1)];        
    end
    
    Cd=zeros(length(r),1);
    distances=pdist(X,distancetype);
    for i = 1:length(r)
        Cd(i)=sum(distances<=r(i))/length(distances);
    end
end