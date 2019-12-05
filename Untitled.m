附录 部分多目标粒子群优化算法源程序
附录1 基本粒子群优化算法源程序
%主程序（最小化问题）
clear all;
t0=clock;

%设置相关参数
c1=1.4962;                                    %学习因子1
c2=1.4962;                                    %学习因子2
w=0.7298;                                     %惯性权重
Tmax=1000;                                    %最大迭代次数
k=10;                                         %决策变量个数
popsize=40;                                   %粒子群中粒子的数目
bounds(1:k,2)=30;                             %设置决策变量取值的上界
bounds(1:k,1)=-30;                            %设置决策变量取值的下界
Global popsize k bounds                       %全局化相关变量

%初始化粒子的位置
range=(bounds(:,2)-bounds(:,1))';    
pop=zeros(popsize,k);
pop(:,1:k)=(ones(popsize,1)*bounds(:,1)')+(ones(popsize,1)*range).*(rand(popsize,numVars));

%初始化粒子的速度
V=randn(popsize,k);

%评价粒子的适应值，并初始化粒子的个体引导者和全局引导者
popeff=evaluation(pop);                     %evaluation为粒子适应值评价函数
Lbest=popeff;                               %Lbest为粒子的个体引导者
[value,index]=sort(Lbest(:,k+1));
Gbest(i,:)=Lbest(index(1),:);

%主循环
GG=[];
for t=1:Tmax
  for i=1:popsize
    for j=1:k  
      V(i,j)=w*V(i,j)+c1*rand*(Lbest(i,j)-pop(i,j))
           +c2*rand*(Gbest(1,j)-pop(i,j));
      pop(i,j)=pop(i,j)+V(i,j);
      if(pop(i,j)<bounds(j,1))
         pop(i,j)=bounds(j,1);
      elseif(pop(i,j)>bounds(j,2))
         pop(i,j)=bounds(j,2);
      end
    end
  end
  popeff=evaluation(pop,k);                    %评价新生粒子的适应值
  for(i=1:1:popsize)
    if popeff(i,k+1)<=Lbest(i,k+1)             %更新粒子的个体引导者
        Lbest(i,:)=popeff(i,:);
end
if popeff(i,k+1)<=Gbest(i,k+1)             %更新粒子的全局引导者
   Gbest(i,:)=popeff(i,:);
end
  end 
  if t/10==ceil(t/10)                           %每隔10代保存一次全局最优点
      GG=[GG;Gbest(k+1)];                     
  end
end

%显示结果
xt=1:10:Tmax
plot(xt,GG(k+1))；                              %展示最优结果的进化曲线
tt1=etime(clock,t0)                             %输出执行时间

% 粒子适应值的计算
function eff=evaluation(pop)
% 函数Rosenbrock，决策变量维数k，取值范围bounds=[-30,30]) 最优值为 F(1)=0 

eff=pop(:,1:k);
x=zeros(1,k);
for i=1:popsize
    x=pop(i,:);d=0;d1=0;
    for j=1:k-1
        d=d+(100*(x(j+1)-x(j)^2)^2+(x(j)-1)^2);
    end
    eff(i,k+1)=d;
    x=zeros(1,k);
end

附录2 第4章多目标骨干粒子群优化算法源程序
%主程序（最小化问题）
clear all;
t0=clock;

%设置相关参数
M=3;                                               %目标函数个数
k=30;                                              %决策变量维数
bounds(1:k,1)=0;                                   %决策变量取值下界
bounds(1:k,2)=1;                                   %决策变量取值上界 
Na=50;                                             %储备集最大容量
popsize=30;                                        %粒子群规模
Tmax=300;                                          %算法最大迭代次数
Global popsize M k bounds Na                       %全局化相关变量

%初始化粒子的位置
range=(bounds(:,2)-bounds(:,1))';    
pop=zeros(popsize,k);
pop(:,1:k)=(ones(popsize,1)*bounds(:,1)')+(ones(popsize,1)*range).*(rand(popsize,numVars));

%初始化粒子的速度
V=randn(popsize,k);

%评价粒子的适应值，并初始化粒子的个体引导者和全局引导者
Lbest=evaluation(pop);                  %初始化粒子的个体引导者
AC=[];                                  %初始化储备集为空集[AC,Gbest]=up_vac(Lbest,AC);            %初始化储备集及全局引导者
 

%主循环
for t=1:Tmax
  pop=up_pop(pop,Gbest,Lbest);          %更新粒子的位置
  EFF=evaluation(pop);                  %评价粒子的适应值
  Lbest=up_Lbest(EFF,Lbest);            %更新粒子的个体引导者
  [AC,Gbest]=up_vac(Lbest,AC);          %更新储备集并选择粒子的全局引导者
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

% 粒子适应值的计算(列举了常用12个测试函数)
function eff=evaluation(pop)

%函数schaffer，变量取值范围[-5,10]，决策变量维数k=1,目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);
  if(x(1)<=1)
    f1=-x(1);
  elseif(x(1)>1 & x(1)<=3)
    f1=x(1)-2;
  elseif(x(1)>3 & x(1)<=4)
    f1=4-x(1);
  else
    f1=x(1)-4;
  end
  f2=(x(1)-5)^2;
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,k);
end

%函数fonseca，变量取值范围[-4,4]，决策变量维数k=8,目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);
  d1=0;d2=0;
  for(p=1:k)
    d1=d1+(x(p)-1/sqrt(k))^2;
    d2=d2+(x(p)+1/sqrt(k))^2;
  end
  f1=1-exp(-d1);
  f2=1-exp(-d2);
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,k);
end

%函数KUR，变量取值范围[-5,5]，决策变量维数k=3,目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);
  f1=-10*exp(-0.2*sqrt(x(1)^2+x(2)^2))-
    10*exp(-0.2*sqrt(x(2)^2+x(3)^2));
  f2=abs(x(1))^0.8+5*sin(x(1)^3)+abs(x(2))^0.8
       +5*sin(x(2)^3)+abs(x(3))^0.8+5*sin(x(3)^3);
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,3);
end

%函数quagliarella，变量取值范围[-5.12,5.12]，决策变量维数k=10，目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);
  d1=0;d2=0;
  for p=1:k
     d1=d1+x(p)^2-10*cos(2*pi*x(p))+10;
     d2=d2+(x(p)-1.5)^2-10*cos(2*pi*(x(p)-1.5))+10;
  end
  f1=sqrt(d1/k);
  f2=sqrt(d2/k);
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,k);
end
 
%函数ZTD1，变量取值范围[0,1]，决策变量维数k=100，目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);d=0;
  for p=2:k
     d=d+x(p);
  end
  g=1+9*d/(k-1);
  f1=x(1);
  f2=g*(1-sqrt(x(1)/g));
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,k);
end
 
%函数ZTD2，变量取值范围[0,1]，决策变量维数k=100，目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);d=0;
  for p=2:k
     d=d+x(p);
  end
  g=1+9*d/(k-1);
  f1=x(1);
  f2=g*(1-(x(1)/g)^2);
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,k);
end

%函数ZTD3，变量取值范围[0,1]，决策变量维数k=100，目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);d=0;
  for j=2:1:size(x,2)
     d=d+x(j);
  end
  g=1+9*d/(k-1);
  f1=x(1);
  f2=g*(1-sqrt(x(1)/g)-(f1/g)*sin(10*pi*f1));
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,k);
end
 
%函数ZTD4，变量取值范围[0,1]，决策变量维数k=30，目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);d=0;
  for j=2:k
     d=d+x(j)^2-10*cos(4*pi*x(j));
  end
  g=1+10*(k-1)+d;
  f1=x(1);
  f2=g*(1-sqrt(x(1)/g));
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,k);
end
 
%函数ZTD6，变量取值范围[0,1]，决策变量维数k=30，目标函数个数M=2
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);d=0;
  for p=2:k
     d=d+x(p);
  end
  g=1+9*(d/(k-1))^(1/4);
  f1=1-exp(-4*x(1))*(sin(6*pi*x(1)))^6;
  f2=g*(1-(f1/g)^2);
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  x=zeros(1,k);
end
 
%函数DTLZ1，变量取值范围[0,1]，决策变量维数k=7，目标函数个数M=3
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);
  g2=0;
  for(j=3:k)
     g2=g2+(x(j)-0.5)^2-cos(20*pi*(x(j)-0.5));
  end
  g=100*(k-2+g2);
  f1=0.5*x(1)*x(2)*(1+g);
  f2=0.5*x(1)*(1-x(2))*(1+g);
  f3=0.5*(1-x(1))*(1+g);
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  eff(i,k+3)=f3;
  x=zeros(1,k);
end
 
%函数DTLZ2，变量取值范围[0,1]，决策变量维数k=12，目标函数个数M=3
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);g=0;
  for(p=3:1:k)
     g=g+(x(p)-0.5)^2;
  end
  f1=(1+g)*cos(pi*x(1)/2)*cos(pi*x(2)/2);
  f2=(1+g)*cos(pi*x(1)/2)*sin(pi*x(2)/2);
  f3=(1+g)*sin(pi*x(1)/2);
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  eff(i,k+3)=f3;
  x=zeros(1,k);
end
 
%函数DTLZ3，变量取值范围[0,1]，决策变量维数k=10，目标函数个数M=3
eff=pop(:,1：k);
x=zeros(1,k);
for i=1:popsize
  x=pop(i,:);g=0;
  for(p=3:1:k)
      g=g+(x(p)-0.5)^2-cos(20*pi*(x(p)-0.5));
  end
  g=100*[k-3+1+g];
  f1=(1+g)*cos(pi*x(1)/2)*cos(pi*x(2)/2);
  f2=(1+g)*cos(pi*x(1)/2)*sin(pi*x(2)/2);
  f3=(1+g)*sin(pi*x(1)/2);
  eff(i,k+1)=f1;
  eff(i,k+2)=f2;
  eff(i,k+3)=f3;
  x=zeros(1,k);
end
 
% pareto支配比较函数
function new_AC=sel_pareto(BEST)
  oldAC=BEST(1,:);
  for(ii=1:size(BEST,1))
      ss1=size(oldAC,1);
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

% 粒子位置更新函数
function pop=up_pop(pop,Gbest,Lbest)
 
pp=exp(-10*t/tt);
for(i=1:1:popsize)
  for(j=1:1:k)
    c1=rand;
    c2=1-c1;
    if(rand<0.5)
      a=(c1*Lbest(i,j)+c2*Gbest(i,j));
      b=abs(Lbest(i,j)-Gbest(i,j));
      pop(i,j)= a + b*randn;
    else
      pop(i,j)=Gbest(i,j);
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

%更新粒子的个体引导者
function  Lbest=up_lbest(EFF,Lbest) 

for(i=1:popsize) 
  cc=0;
  for(j=1:M)
    if(EFF(i,k+j)<=Lbest(i,k+j))
      cc=cc+1;
    end
  end
  if(cc~=0)
    Lbest(i,:)=EFF(i,:);               
  end
end

%更新算法外部储备集，以及粒子的全局引导者
function [newAC,gbest]=up_vac(AC,oldAC)  

For i=1:size(AC,1)
    oldAC=up_vac0(AC(i,:),oldAC);  
end
newAC=oldAC;
crowd_value=calcul_crowd(newAC);
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

%更新储备集中的元素
function new_AC=up_vac0(par_eff,oldAC)  
ss1=size(oldAC,1);
For i=1:ss1
   bb1=0;bb2=0;
   For j=1:M
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
 if(ss2>Na)
    for(i=1:M)
         LIM_f(i,2)=max(new_AC(:,k+i));
         LIM_f(i,1)=min(new_AC(:,k+i));
    end
    DD=[];
    deep=[];
    for(i=1:M)
        [val,ind]=sort(new_AC(:,k+i));
        for(j=1:Na+1)
            if(j==1 | j==Na+1 )
                DD(ind(j),i)=inf;
            else                
                DD(ind(j),i)=(new_AC(ind(j+1),k+i)
                      -new_AC(ind(j-1),k+i))/(LIM_f(i,2)-LIM_f(i,1));  
            end
        end
    end
    for(jj=1:Na+1)
       deep(jj)=sum(DD(jj,:));
    end
    [val,ind]=sort(deep);
     new_AC(ind(1),:)=[];
 end       

%计算储备集中元素的拥挤距离值
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

