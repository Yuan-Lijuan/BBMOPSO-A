 function  LBEST=ppso1(EFF,LBEST,popsize,M,k)

%����Ⱥ�����Ӹ��弫ֵ�ĸ��£�������������������������������


   
  for(i=1:popsize) 
      cc=0;
      for(j=1:M)
        if(EFF(i,k+j)<=LBEST(i,k+j))
            cc=cc+1;
        end
      end
      if(cc~=0)
          LBEST(i,:)=EFF(i,:);               
      end

  end