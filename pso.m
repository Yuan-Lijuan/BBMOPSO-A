
% function  [AC,comp_solut,best_cost,best_emiss]=pso
%    clear

%{
张老师的源代码

M=2;
k=8;
bounds(1:k,1)=0;

bounds(1:k,2)=1;


arch_size=50;


t0=clock;
popsize=30;


tt=300;


pop=initialize(popsize,bounds,k);
LBEST=evaluation(pop,k);
AC=sel_pareto(LBEST,M,k);
[AC,gbest]=up_vac(AC,AC,arch_size,popsize,M,k);

t=0;
while(t<=tt)
    pop=ppso(pop,gbest,LBEST,M,k,popsize,t,tt,bounds);
    EFF=evaluation(pop,k);
    LBEST=ppso1(EFF,LBEST,popsize,M,k);
    [AC,gbest]=up_vac(LBEST,AC,arch_size,popsize,M,k);
    t=t+1
end
if(M==3)
    figure;
    plot3(AC(:,k+1),AC(:,k+2),AC(:,k+3),'black.')
else
    figure;
    plot(AC(:,k+1),AC(:,k+2),'black.')
    hold on

end
hold on
ttt=etime(clock,t0)
%}
clear
clc
%MyBBMOPSO
M=2;
k=12;
bounds(1:k,1)=[0,0,0,2,0.0001,0.0001,0.1,0.1,6,10,18,24];
%将决策变量范围转换模型到空间中的范围
bounds(1:k,2)=[360,3.6,3.9,6,0.7,0.1,1,1,12,18,23,28];
%boundu=[360,3.6,3.9,1,1,1,12,18,23,28];
%boundl=[0,0,0,0.1,0.1,0.1,6,10,18,24];
arch_size=50;
t0=clock;
popsize=50;
tt=20;
pop=initialize(popsize,bounds,k)


LBEST=evaluation(pop,k)

 AC=sel_pareto(LBEST,M,k);
 [AC,gbest]=up_vac(AC,AC,arch_size,popsize,M,k);

 t=0;
 while(t<=tt)
%   pop=ppso(pop,gbest,LBEST,M,k,popsize,t,tt,bounds); 
    pop=ppso(pop,gbest,LBEST,M,k,popsize,t,tt,bounds,AC);
    EFF=evaluation(pop,k);                             
    LBEST=ppso1(EFF,LBEST,popsize,M,k);                
    [AC,gbest]=up_vac(LBEST,AC,arch_size,popsize,M,k);
    t=t+1
end
if(M==3)
    figure;
   plot3(AC(:,k+1),AC(:,k+2),AC(:,k+3),'black.');
else
  figure;
  plot(AC(:,k+1),AC(:,k+2),'black.')
  hold on
%
 end
 hold on
 ttt=etime(clock,t0)
% 
% 
% 
% 
% 
% 







