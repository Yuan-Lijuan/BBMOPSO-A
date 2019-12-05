 function new_AC=sel_pareto(BEST,M,k)
%用来进行储备集更新
  % psize            the number of populations in population
  %tmobsize           he dimension of multiobjective 
  %OBM           the values matrix of the population in the space of multiobjective
  %PAM              储备集中个体对应的目标空间矩阵  
  oldAC=BEST(1,:);
  for(ii=1:size(BEST,1))   %返回BEST的行数
      ss1=size(oldAC,1);   %返回old的行数
      for(i=1:ss1)
          bb1=0;bb2=0;
          for(j=1:M)
              aa1=oldAC(i,k+j);
              aa2=BEST(ii,k+j);
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
              BEST(ii,:)=inf;
              break;
          elseif(bb2~=0 & bb1==0)
              BEST(ii,:)=inf;
              break;
          end
      end
      sum_AC=[BEST(ii,:);oldAC];
      new_AC=[];
      for(i=1:size(sum_AC,1))
          if(sum_AC(i,1)~=inf)
              new_AC=[new_AC;sum_AC(i,:)];
          end
      end
      oldAC=new_AC;
  end
  
  