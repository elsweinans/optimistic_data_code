%% eg calc_corrdims("int_resp","b",5,1,[0.05 0.5], 6, 'euclidean',51840)


function [slopes] = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n)
    
    ND = load(group + "_group");
    try
    inds = ND.(group + "_group"); % credits to Victor Azizi
    catch ERROR
      disp("field name in file should match filename")
      rethrow(ERROR)
    end
 
    ks=length(inds);
    slopes=zeros(ks,1);
    
    for k=1:ks
       groupind=inds(k);
       if group == "healthy"
           filename=group + "/" + groupind + ".csv";
       else
           filename=group + "/" + groupind + "_" + timing + ".csv";
       end
       data = readtable(filename);
       data1=data.x; 
       idxact=find(data1>0.05);
       data2=data1(idxact(1):idxact(end));
       
       mvavg=smoothdata(data2,'gaussian',smoothwindow);
              
       L=floor(length(mvavg)/avgwindow);
       mvavg2=mean(reshape(mvavg(1:L*avgwindow),[avgwindow,L]),1)';
       
       if n<length(mvavg2)
           disp(filename)
           mvavg2=mvavg2(1:n); 
       end
       
 
       
%        figure
%        plot(mvavg2)
       
       nrsteps=10;
       rs=logspace(log(ranger(1)*std(mvavg2))/log(10),log(ranger(2)*std(mvavg2))/log(10),nrsteps); %mvavg2 ipv data1
       Cds=calc_Cd(mvavg2,D,distancetype,rs);     
       slopefit=polyfit(log(rs'),log(Cds),1);
       slopes(k)=slopefit(1);
    end
end


