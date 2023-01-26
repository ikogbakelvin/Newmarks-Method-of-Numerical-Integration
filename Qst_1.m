clc
clear all
close all

method='Newmarks_method';

%%Using Newmarkmethod (Constant and Linear acceleration)
%%Example 7.1 page 191.Response of damped system using Newmarkmethod.

%%Earthquake groundmotion
ug=xlsread('JMA Kobe EW.xlsx','B2:B7578');
ug=ug.';

%%Given parameters
m=1; %%assume
Time=1;
dt=0.02;
Time_history=0:dt:40.96;
wn=10;
zeta=0.02;%Damping ratio
k=wn*wn*m;
wd=wn*sqrt(1-zeta^2);
c=2*zeta*m*wn; 
N=length(Time_history);

betav= [1/4 1/6];
for i=1:2
    beta = betav(i);
    gamma = 1/2;
end


%%Initial Calculation

a1=(1/(beta*(dt^2)))*m+(gamma/(beta*dt))*c;
a2=(1/(beta*dt))*m+((gamma/beta)-1)*c;
a3=((1/(2*beta))-1)*m+dt*((gamma/(2*beta))-1)*c;
kcap=k+a1;

%%Initial Calculation
x0=0;
xd0=0;
xdd0=(-m*ug-c*xd0-k*x0)/m;

x(1)=x0;
xd(1)=xd0;
xdd(1)=-ug(1);

%%Repetition for ith time steps


for i=1:N
    t(i+1)=dt*(i-1);
    pcap(i+1)=-m*ug(i+1)+a1*x(i)+a2*xd(i)+a3*xdd(i);
    x(i+1)=pcap(i+1)/kcap;
    xd(i+1)=(gamma/(beta*dt))*(x(i+1)-x(i))+(1-(gamma/beta))*xd(i)+dt*(1-(gamma/(2*beta)))*xdd(i);
    xdd(i+1)=(1/(beta*(dt^2)))*(x(i+1)-x(i))-(1/(beta*dt))*xd(i)-((1/(2*beta))-1)*xdd(i);
end

for i=1:N+1
    xddt(i)=xdd(i)+ ug(i);
    
end


t=t.';
x=x.';
xd=xd.';
xdd=xdd.';



filename='JMA Kobe EW';
v=[t,x,xd,xdd];
str={'time','disp','vel','acc'};
x_c=num2cell(v);
data=[str;x_c];
xlswrite([filename,'.xlsx'],data,method);

subplot(2,1,1)
plot(t, x)
grid on
xlabel('time (sec)')
ylabel('Displacement (m)')
title('Displacement Response')

subplot(2,1,2)
plot(t, xd)
grid on
xlabel('time (sec)')
ylabel('Velocity (m/s)')
title('Velocity Response')

figure
plot(t, xdd)
grid on
xlabel('time (sec)')
ylabel('Acceleration (m/s2)')
title('Acceleration History')

subplot(2,1,2)
plot(t, xddt)
grid on
xlabel('time (sec)')
ylabel('Acceleration (m/s2)')
title('Total Acceleration Response')


