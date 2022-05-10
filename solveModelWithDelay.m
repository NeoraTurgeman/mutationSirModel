function solveModelWithDelay()
global t;
t=linspace(0,200,1);
lags=[14 25];
global n;
n= 9000000;

beta_empty_1 = 0.4  ;
	beta_empty_2 = 0.3 ;
	beta_1_2 =     0.25  ;
	gama_empty_1 = 0.1 ;
	gama_empty_2 = 0.1 ;
	gama_1_2 = 0.125  ;
	rho_empty_1 = 0.98 ;
	rho_empty_2 = 0.96 ;
	rho_1_2 = 0.99 ;

tspan = [0,200];
sol = dde23(@ddefun, lags, @history , tspan);
 yint = deval(sol,tspan);
%yint;
max=154;
for i = 1:154
if  sol.y(1,i)<0
    for j = i:max
    sol.y(1,i)=0;
    s2=sol.y(2,i);
    sol.y(2,i)=sol.y(2-1,i)-gama_empty_1*sol.y(2,i-1);
    sol.y(3,i)=sol.y(3-1,i)-gama_empty_2*sol.y(3,i-1);
    sol.y(4,i)=sol.y(4,i)+gama_1_2*(-s2+sol.y(4,i-1));
    end
end
if  sol.y(2,i)<0
    for j = i:max
    s2=sol.y(2,i);
    sol.y(2,i)=0;
    sol.y(4,i)=sol.y(4,i)+gama_1_2*(-s2+sol.y(4,i-1));
    end
end
if  sol.y(4,i)<0
    for j = i:max
        sol.y(4,i)=0;
        s6= sol.y(6,i);
    sol.y(6,i)=sol.y(6-1,1)-gama_1_2*sol.y(6,i-1);
    sol.y(7,i)=gama_1_2*rho_1_2*sol.y(6,i);
    sol.y(8,i)=sol.y(8,i)+gama_1_2*(1-rho_1_2)*(-s6+sol.y(6,i));
        end
end
end


%y=deval(sol,t);
figure;
plot(sol.x,sol.y,'-');
xlabel('Time t');
ylabel('Solution y');
legend('R_0','R_0I_1','R_0I_2', 'R_1', 'R_2', 'R_1I_2', 'R_{1,2}', 'D', 'Location','NorthWest');

end

function dydt = ddefun(t,y,Z)
global n;
n= 9000000;
    beta_empty_1 = 0.3  ;
	beta_empty_2 = 0.25 ;
	beta_1_2 =     0.25  ;
	gama_empty_1 = 0.1 ;
	gama_empty_2 = 0.1 ;
	gama_1_2 = 0.125  ;
	rho_empty_1 = 0.98 ;
	rho_empty_2 = 0.96 ;
	rho_1_2 = 0.99 ;
    lags=[14 25];
    ylag1 = Z(:,1);
    ylag2 = Z(:,2);

  dydt = [%R0
          -y(1)*(beta_empty_1 / n)*y(2)-ylag1(1)*(beta_empty_2 / n)*(ylag1(3)+ylag1(6));
          %R0I1
          (beta_empty_1 / n)*y(2)*y(1)-gama_empty_1*y(2);
          %R0I2
          (beta_empty_2 / n)*(ylag1(3)+ylag1(6))*ylag1(1)-(gama_empty_2*ylag1(3));
          %R1
          gama_empty_1*rho_empty_1*y(2)-((beta_1_2 / n)*(ylag1(3)+ylag1(6))*ylag1(4));
          %R2
          gama_empty_2*rho_empty_2*ylag1(3);
          %R1I2
          (beta_1_2 / n)*(ylag1(3)+ylag1(6))*ylag1(4)-gama_1_2*ylag1(6);
          %R12
          gama_1_2*rho_1_2*ylag1(6);
          %D
          gama_empty_1*(1-rho_empty_1)*y(2)+gama_empty_2*(1-rho_empty_2)*ylag1(3)+gama_1_2*(1-rho_1_2)*ylag1(6);]
end

function s = history(t)
  s = [9000000;1300;1300;0;0;0;0;54];
end