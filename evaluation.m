
function eff=evaluation(pop,k)


% %______________________KUR-------------------------(-5,5)
% % % 
% 
% for(i=1:1:size(pop,1))
%     for(j=1:1:size(pop,2))
%  eff(i,j)=pop(i,j);
% end
% end
% x=zeros(1,3);
% for( i=1:1:size(pop,1))
%     x=pop(i,:);
% f1=-10*exp(-0.2*sqrt(x(1)^2+x(2)^2))-10*exp(-0.2*sqrt(x(2)^2+x(3)^2));
% f2=abs(x(1))^0.8+5*sin(x(1)^3)+abs(x(2))^0.8+5*sin(x(2)^3)+abs(x(3))^0.8+5*sin(x(3)^3);
% 
%     eff(i,4)=f1;
%     eff(i,5)=f2;
%                                                                                                                                                                                                                                         
%     x=zeros(1,3);
% end

%%%%%%%%%%%%%%----------function quagliarella---2  ---[-5.12,5.12]---
% 
% for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                   d1=0;d2=0;
%                   for(p=1:1:k)
%                      d1=d1+x(p)^2-10*cos(2*pi*x(p))+10;
%                      d2=d2+(x(p)-1.5)^2-10*cos(2*pi*(x(p)-1.5))+10;
%                   end
%                   
%                   f1=sqrt(d1/k);
%                   f2=sqrt(d2/k);
%                    eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end

% 
% 
% %--------ZTD1------------------------------------------
% % %   
for(i=1:1:size(pop,1))
    for(j=1:1:size(pop,2))
        eff(i,j)=pop(i,j);
    end
end
%x=zeros(1,k);
for( i=1:1:size(pop,1))
%     x=pop(i,:);d=0;
%     for(p=2:1:k)
%         d=d+x(p);
%     end
%     g=1+9*d/(k-1);


   % x=[0 1 2];
    x=pop(i,:);
    printf=x;
    dlmwrite('C:\EnergyPlusV8-1-0\ExampleFiles\test1\TwelveVariables.txt',x,' ');
%     dlmwrite('FiveVarious.txt',x,' ')
    cmd = 'C:\EnergyPlusV8-1-0\ExampleFiles\test1\6\Debug\6.exe';
    system(cmd);
    system('C:\EnergyPlusV8-1-0\RunEPlus 11yue22 CHN_Jiangsu.Xuzhou.580270_CSWD');
    f1=xlsread('C:\EnergyPlusV8-1-0\ExampleFiles\test1\11yue22Table.csv','A15:C15');
    f2=xlsread('C:\EnergyPlusV8-1-0\ExampleFiles\test1\11yue22Table.csv','A167:C167');
    a=size(f2);
    if (a==0)
        f2=xlsread('C:\EnergyPlusV8-1-0\ExampleFiles\test1\11yue22Table.csv','A166:C166');
   else
        f2=xlsread('C:\EnergyPlusV8-1-0\ExampleFiles\test1\11yue22Table.csv','A167:C167');
    end
    
    eff(i,k+1)=f1;
    eff(i,k+2)=f2;
    
%     x=zeros(1,k);
end





% 
% 
% % %  %____________________ZDT2_____________________________
%  for(i=1:1:size(pop,1))
%                  for(j=1:1:size(pop,2))
%                   eff(i,j)=pop(i,j);
%                         end
%            end 
%           x=zeros(1,k);
%           for( i=1:1:size(pop,1))
%                        x=pop(i,:);d=0;
%                   for(p=2:1:k)
%                       d=d+x(p);
%                    end
%                    g=1+9*d/(k-1);
%                      f1=x(1);
%                    f2=g*(1-(x(1)/g)^2);
%                       eff(i,k+1)=f1;
%                       eff(i,k+2)=f2;
%                    x=zeros(1,k);
%        end
% 


% 
  %----------------------------ZDT3__________________________     
% % %  
% % 
% for(i=1:1:size(pop,1))
%     for(j=1:1:size(pop,2))
%  eff(i,j)=pop(i,j);
% end
% end
% x=zeros(1,k);
% for( i=1:1:size(pop,1))
%     x=pop(i,:);
%     d=0;
%     for(j=2:1:size(x,2))
%         d=d+x(j);
%     end
%     g=1+9*d/(k-1);
%    f1=x(1);
%    f2=g*(1-sqrt(x(1)/g)-(f1/g)*sin(10*pi*f1));
%     eff(i,k+1)=f1;
%     eff(i,k+2)=f2;
%     x=zeros(1,k);
% end
% 
%  %----------------------------ZDT4__________________________     
% % % %  
% % 
% for(i=1:1:size(pop,1))
%     for(j=1:1:size(pop,2))
%  eff(i,j)=pop(i,j);
% end
% end
% x=zeros(1,k);
% for( i=1:1:size(pop,1))
%     x=pop(i,:);
%     d=0;
%     for(j=2:1:k)
%         d=d+x(j)^2-10*cos(4*pi*x(j));
%     end
%     g=1+10*(k-1)+d;
%    f1=x(1);
%    f2=g*(1-sqrt(x(1)/g));
%     eff(i,k+1)=f1;
%     eff(i,k+2)=f2;
%     x=zeros(1,k);
% end




% % %----------------------------ZDT6______________________
% 
% % % 
% for(i=1:1:size(pop,1))
%     for(j=1:1:size(pop,2))
%  eff(i,j)=pop(i,j);
% end
% end
% x=zeros(1,k);
% for( i=1:1:size(pop,1))
%     x=pop(i,:);d=0;
%     for(p=2:1:k)
%         d=d+x(p);
%     end
%     g=1+9*(d/(k-1))^(1/4);
%    f1=1-exp(-4*x(1))*(sin(6*pi*x(1)))^6;
%    f2=g*(1-(f1/g)^2);
%     eff(i,k+1)=f1;
%     eff(i,k+2)=f2;
%     x=zeros(1,k);
% end

% 
% % %---------------DTLZ1(3 wei)-----------(0,1)
% for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                      g2=0;
%                      for(j=3:k)
%                          g2=g2+(x(j)-0.5)^2-cos(20*pi*(x(j)-0.5));
%                      end
%                      g=100*(k-2+g2);
%                      f1=0.5*x(1)*x(2)*(1+g);
%                   f2=0.5*x(1)*(1-x(2))*(1+g);
%                   f3=0.5*(1-x(1))*(1+g);
%                  eff(i,k+1)=f1;
%                  eff(i,k+2)=f2;
%                  eff(i,k+3)=f3;               
%                    x=zeros(1,k);
%            end
% %            
%            

% % %---------------DTLZ2(3 wei)-----------(0,1)
% 
% 
% 
% 
% 
% for(i=1:1:size(pop,1))
%     for(j=1:1:size(pop,2))
%  eff(i,j)=pop(i,j);
% end
% end
% x=zeros(1,k);
%   for( i=1:1:size(pop,1))
%     x=pop(i,:);g=0;
%     for(p=3:1:k)
%         g=g+(x(p)-0.5)^2;
%     end
% f1=(1+g)*cos(pi*x(1)/2)*cos(pi*x(2)/2);
% f2=(1+g)*cos(pi*x(1)/2)*sin(pi*x(2)/2);
%  f3=(1+g)*sin(pi*x(1)/2);   
%     eff(i,k+1)=f1;
%     eff(i,k+2)=f2;
%     eff(i,k+3)=f3;
%     x=zeros(1,k);
% end
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-ZDT3-----------------

% 
% for(i=1:1:size(pop,1))
%     for(j=1:1:size(pop,2))
%  eff(i,j)=pop(i,j);
% end
% end
% x=zeros(1,k);
%   for( i=1:1:size(pop,1))
%     x=pop(i,:);g=0;
%     for(p=3:1:k)
%         g=g+(x(p)-0.5)^2-cos(20*pi*(x(p)-0.5));
%     end
%     g=100*[k-3+1+g];
% f1=(1+g)*cos(pi*x(1)/2)*cos(pi*x(2)/2);
% f2=(1+g)*cos(pi*x(1)/2)*sin(pi*x(2)/2);
%  f3=(1+g)*sin(pi*x(1)/2);   
%     eff(i,k+1)=f1;
%     eff(i,k+2)=f2;
%     eff(i,k+3)=f3;
%     x=zeros(1,k);
% end
% end





% 

%--------fonseca(2) -------------[-4,4]-----------------------------
%   
%    for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                   d1=0;
%                   d2=0;
%                   for(p=1:k)
%                      d1=d1+(x(p)-1/sqrt(k))^2;
%                      d2=d2+(x(p)+1/sqrt(k))^2;
%                   end
%                   f1=1-exp(-d1);
%                   f2=1-exp(-d2);
%                  eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end



% 
% % ----------function MOP6---2  ---[0,1]---
% a=2;
% q=4;
% for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                   f1=x(1);
%                   f2=(1+10*x(2))*(1-(x(1)/(1+10*x(2)))^a-x(1)/(1+10*x(2))*sin(2*pi*q*x(1)));
%                    eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end
% 















% % % ----------function Binh(1)---2  ---[-5,10]---
%    for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                   f1=x(1)^2+x(2)^2;
%                   f2=(x(1)-5)^2+(x(2)-5)^2;
%                    eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end






%--------fonseca -------------none-----------------------------
%   
%    for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                      f1=1-exp(-(x(1)-1)^2-(x(2)+1)^2);
%                   f2=1-exp(-(x(1)+1)^2-(x(2)-1)^2);
%                  eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end



%----------function LIS---2  ---[-5,10]---
%    for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                   f1=(x(1)^2+x(2)^2)^(1/8);
%                   f2=((x(1)-0.5)^2+(x(2)-0.5)^2)^(1/4);
%                    eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end
% 





%----------function POLONI---2  ---[-pi,pi]---
%  A1=0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
%  A2=1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);
% for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                    B1=0.5*sin(x(1))-2*cos(x(1))+sin(x(2))-1.5*cos(x(2));
%                    B2=1.5*sin(x(1))-cos(x(2))+2*sin(x(2))-0.5*cos(x(2));
%                   
%                   f1=-(1+(A1-B1)^2+(A2-B2)^2);
%                   f2=-((x(1)+3)^2+(x(2)+1)^2);
%                    eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end








% % % ----------function schaffer---2  ---[-5,10]---
% 
% for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                   if(x(1)<=1)
%                     f1=-x(1);
%                   elseif(x(1)>1 & x(1)<=3)
%                       f1=x(1)-2;
%                   elseif(x(1)>3 & x(1)<=4)
%                       f1=4-x(1);
%                   else
%                       f1=x(1)-4;
%                   end
%                   f2=(x(1)-5)^2;
%                    eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end
% 

% % % % ----------function vicini---3  ---s [-2,2]---
%    for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                    eff(i,k+1)=x(1)^2+(x(2)-1)^2;
%                    eff(i,k+2)=x(1)^2+(x(2)+1)^2+1;
%                    eff(i,k+3)=x(2)^2+(x(1)-1)^2+2;
%                    x=zeros(1,k);
%            end
% 


% % % % % ----------function vicini---3  ---s [-4,4]---
%    for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                    eff(i,k+1)=((x(1)-2)^2)/2+((x(2)+1)^2)/13+3;
%                    eff(i,k+2)=((x(1)+x(2)-3)^2)/36+((x(2)-x(1)+2)^2)/8-17;
%                    eff(i,k+3)=(x(1)-2*x(2)-1)^2/175+((2*x(2)-x(1))^2)/17-13;
%                    x=zeros(1,k);
%            end



% % % % % ----------function vicini---3  ---s [-3,3]---
%    for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                    eff(i,k+1)=0.5*(x(1)^2+x(2)^2)+sin(x(1)^2+x(2)^2);
%                    eff(i,k+2)=(3*x(1)-2*x(2)+4)^2/8+(x(1)-x(2)+1)^2/27+15;
%                    eff(i,k+3)=1/(x(1)^2+x(2)^2+1)-1.1*exp(-x(1)^2-x(2)^2);
%                    x=zeros(1,k);
%            end
% 








%_____________________________________Bimodal problem proposed by Deb [1]: 

% for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                    g=2.0-exp(-((x(2)-0.2)^2)/0.004)-0.8*exp(-((x(2)-0.6)^2)/0.4);
%                 
%                      f1=x(1);
%                   f2=g/x(1);
%                  eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    x=zeros(1,k);
%            end


%------binh(3)-----3 wei---[10^-6,10^6]

% 
%    for(i=1:1:size(pop,1))
%        for(j=1:1:size(pop,2))
%         eff(i,j)=pop(i,j);
%           end
%        end
%                x=zeros(1,k);
%             for( i=1:1:size(pop,1))
%                   x=pop(i,:);
%                   f1=x(1)-10^6;
%                   f2=x(2)-2*10^(-6);
%                   f3=x(1)*x(2)-2;
%                    eff(i,k+1)=f1;
%                    eff(i,k+2)=f2;
%                    eff(i,k+3)=f3;
%                    x=zeros(1,k);
%            end
% 
