
function pop=ppso(pop,gbest,LBEST,M,k,popsize,t,tt,bounds,AC)
 A=AC(:,13);
 B=AC(:,14);
 f1_max=max(A);
 f2_max=max(B);
 f1_min=min(A);
 f2_min=min(B);
 x_low=bounds(:,1);
 x_up=bounds(:,2);
 
pp=exp(-10*t/tt);
for i=1:1:popsize
    f1_gbest= gbest(i,13);
    f1_LBEST= LBEST(i,13);
    f2_gbest=gbest(i,14);
    f2_LBEST=LBEST(i,14);
    d=f1_LBEST-f1_gbest;
    e=f1_max-f1_min;
    f=f2_LBEST-f2_gbest;
    g=f2_max-f2_min;
    h=abs(d/e)+abs(f/g);
    pro_d=0.5*(1-(1/M)*h);
    for j=1:1:k 
        c=rand;
%         x_up=x_up(j);
%         x_low=x_low(j);
      if(pro_d<c)
         Q(j)=0;
      else
        Q(j)=(x_up(j)-x_low(j))*exp(-5*t/tt);
      end
        c1=rand;
        c2=1-c1;
        if(rand<0.5)
            a=(c1*LBEST(i,j)+c2*gbest(i,j));
%             b=abs(LBEST(i,j)-gbest(i,j));
         b=abs(LBEST(i,j)-gbest(i,j))+Q(j);
            pop(i,j)= a + b*randn;
        else
            pop(i,j)=gbest(i,j);
        end
        if(pop(i,j)<bounds(j,1))
            pop(i,j)=bounds(j,1);
        elseif(pop(i,j)>bounds(j,2))
            pop(i,j)=bounds(j,2);
        end
    end

    if(pp>rand)
        aa=ceil(k*rand);
        rang=(bounds(aa,2)-bounds(aa,1))*pp;
        pop(i,aa)=pop(i,aa)+randn*rang;
        if(pop(i,aa)<bounds(aa,1))
            pop(i,aa)=bounds(aa,1);
        elseif(pop(i,aa)>bounds(aa,2))
            pop(i,aa)=bounds(aa,2);
        end
    end

end