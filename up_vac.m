function [newAC,gbest]=up_vac(AC,oldAC,arch_size,popsize,M,k)  
%------

for(i=1:size(AC,1))
    oldAC=up_vac0(AC(i,:),oldAC,M,arch_size,k);  
end
newAC=oldAC;
crowd_value=calcul_crowd(newAC,M,k);

gbest_set=newAC;
gbest_crowd_val=crowd_value;
g_size=size(gbest_set,1);
while g_size>popsize
    [val,ind]=sort(gbest_crowd_val);
    gbest_set(ind(1),:)=[];
    gbest_crowd_val(ind(1))=[];
    g_size=size(gbest_set,1);
end
for(i=1:popsize)
    a1=ceil(g_size*rand);
    a2=ceil(g_size*rand);
    if(gbest_crowd_val(a1)>=gbest_crowd_val(a2))
        gbest(i,:)=gbest_set(a1,:);
    else
        gbest(i,:)=gbest_set(a2,:);
    end
end





function new_AC=up_vac0(par_eff,oldAC,M,arch_size,k)  

ss1=size(oldAC,1);
for(i=1:ss1)
   bb1=0;bb2=0;
   for(j=1:M)
      aa1=oldAC(i,k+j);
      aa2=par_eff(1,k+j);
      if(aa2<aa1)
          bb1=bb1+1;
      elseif(aa2==aa1)
          bb2=bb2+1;
      end
   end
      if(bb1==M)
         oldAC(i,:)=inf;
      elseif(bb2>0 & bb1==M-bb2)
          oldAC(i,:)=inf;
      elseif(bb1+bb2==0)
          par_eff(1,:)=inf;
          break;
      elseif(bb2~=0 & bb1==0)
         par_eff(1,:)=inf;
         break;
      end    
  end

 sum_AC=[par_eff;oldAC];
 new_AC=[];
 for(i=1:ss1+1)
     if(sum_AC(i,1)~=inf)
         new_AC=[new_AC;sum_AC(i,:)];
     end
 end
 
 
ss2=size(new_AC,1);
 if(ss2>arch_size)
    for(i=1:M)
         LIM_f(i,2)=max(new_AC(:,k+i));
         LIM_f(i,1)=min(new_AC(:,k+i));
    end
    DD=[];
    deep=[];
    for(i=1:M)
        [val,ind]=sort(new_AC(:,k+i));
        for(j=1:arch_size+1)
            if(j==1 | j==arch_size+1 )
                DD(ind(j),i)=inf;
            else                
               DD(ind(j),i)=(new_AC(ind(j+1),k+i)-new_AC(ind(j-1),k+i))/(LIM_f(i,2)-LIM_f(i,1));  
            end
        end
    end
    for(jj=1:arch_size+1)
       deep(jj)=sum(DD(jj,:));
    end
    [val,ind]=sort(deep);
     new_AC(ind(1),:)=[];
 end       

 
 
 
% function new_AC=up_vac0(new_par,oldAC,M,arch_size,k)
% ee=0.05;
% 
% ss1=size(oldAC,1);
% 
% B_oldAC=[];
% for(i=1:ss1)
%     for(j=1:M)
%         B_oldAC(i,j)=floor(log10(oldAC(i,k+j))/log10(1+ee));
%     end
% end
% for(j=1:M)
%        B_par(j)=floor(log10(new_par(k+j))/log10(1+ee));
% end
% for(i=1:ss1)
%    bb1=0;bb2=0;
%    for(j=1:M)
%       aa1=B_oldAC(i,j);
%       aa2=B_par(1,j);
%       if(aa2<aa1)
%           bb1=bb1+1;
%       elseif(aa2==aa1)
%           bb2=bb2+1;
%       end
%    end
%       if(bb1==M)
%          oldAC(i,:)=inf;
%       elseif(bb2>0 & bb1==M-bb2)
%           oldAC(i,:)=inf;
%       elseif(bb1+bb2==0)
%           new_par(1,:)=inf;
%           break;
%       elseif(bb2~=0 & bb1==0)
%           new_par(1,:)=inf;
%          break;
%       end    
%   end
%   
%  sum_AC=[new_par;oldAC];
%  new_AC=[];
%  for(i=1:ss1+1)
%      if(sum_AC(i,1)~=inf)
%          new_AC=[new_AC;sum_AC(i,:)];
%      end
%  end
%  
%  
 
 
function crowd_value=calcul_crowd(new_AC,M,k)


s2=size(new_AC,1);
if(s2>=2)
for(i=1:M)
      LIM_f(i,2)=max(new_AC(:,k+i));
      LIM_f(i,1)=min(new_AC(:,k+i));
end 
DD=[];
crowd_value=[];
    for(i=1:M)
        [val,ind]=sort(new_AC(:,k+i));
        for(j=1:s2)
            if(j==1 )
                DD(ind(j),i)=inf;
            elseif(j==s2)
                DD(ind(j),i)=inf;
            else                
                DD(ind(j),i)=(new_AC(ind(j+1),k+i)-new_AC(ind(j-1),k+i))/(2*(LIM_f(i,2)-LIM_f(i,1)));  
            end
        end
    end
    for(jj=1:s2)
       crowd_value(jj)=sum(DD(jj,:));
    end
else
    crowd_value(1)=1;
end
 
