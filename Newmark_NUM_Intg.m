clc
clear all
close all

method='Newmarks_method';

%%Question 7
%%Using Newmarkmethod (Constant and Linear acceleration)
%%Given parameters

TotalTime=1;
dt=0.02;
Tt=0:dt:TotalTime;
N=length(Tt);

%%Define the Degree of freedom
df=input('Input the degree of freedom\n');

%%Damping ratio for the first and second mode
zeta=[0.05;0.05;0;0];     zeta1=[0.05;0.05];

%%Define the Mass matrix
m=[5 0 0 0;0 5 0 0; 0 0 5 0; 0 0 0 10];

[n,ns]=size(m);

%%Define the Stiffness matrix

k=[400 -200 0 0;-200 400 -200 0;0 -200 400 -200;0 0 -200 600];
   

fprintf('The Modal Mass matrix and modal stiffness matrix of structure are as follow \n')
%%Determine the mode shape and Natural freqency matrix
[V,d]=eig(k,m);

fprintf('The Natural frequencies of the structure are as follow \n')
wn=sqrt(diag(d));
disp(wn);

fprintf('The modal matrix of the structure \n')
disp(V);

%%Determine the modal Damping matrix; C=damping matrix
%%M=modal mass matrix; K=modal stiffnees matrix
Mn=V'*m*V;
disp(Mn);
M1=inv(Mn)*m*V;
M2=inv(Mn)*V'*m;
fprintf('The Modal Stiffness matrix of structure are as follow \n')
Kn=V'*k*V;
disp(Kn);
fprintf('The Modal Damping matrix of structure are as follow \n')
Cn=2*zeta.*Mn*wn;
disp(Cn)
fprintf('The Damping matrix of structure are as follow \n')
C=M1*Cn.*M2;
disp(C)

%%Using Rayleigh Damping Matrix
%%Natural frequencies
W=[1/wn(1) wn(1);
    1/wn(2) wn(2)];
deltaW=det(W);
adjW=adjoint(W);
%%Determine the mass and stiffness proportion a0 and a1
A=(1/deltaW)*adjW*2*zeta1;

fprintf('The mass proportion a0 of Rayleigh damping is\n')
a0=A(1);
disp(a0);
fprintf('The stiffness proportion a1 of Rayleigh damping is\n')
a_1=A(2);
disp(a_1);

%%Determine the modal Damping matrix; C=damping matrix
fprintf('The Daping matrix for Rayleigh damping is\n')
c=a0*m+a_1*k;
disp(c)


%%Determine the Natural Time Period of the structure
for i=1:n
   
    T(i)=2*pi/wn(i);
     
%%Normalizing mode shape vectors
    V(:,i)=V(:,i)/V(n,i);
   
end

%%Natural Period of the structure (sec)
fprintf('Natural Time Period of structure (sec) \n')
disp(T);
%%Modal Shape Matrix
fprintf('The Normalized Modal matrix of structure are as follow \n')
disp(V);

%%Ploting of Modes Shape
%%Assume Height of the Structure

Height=[0;3;6;9;12];

for i=1:n
   
    subplot(1,n,i)
    plot([0; V(:,i)],Height);
    ylabel('Height of the structure','Fontsize',12);
    title(['Mode Shape ',num2str(i)],'Fontsize',12);
   
end

%%Earthquake groundmotion
Time=40.96;
dt=0.02;
t=0:dt:Time;
Z=length(t);
ug=xlsread('JMA Kobe EW.xlsx','B3:B2051');

uddg=ug;

for i=1:n
   
    f(:,i)=-uddg*m(i,i);
   
    for i=1
   
    x(:,i)=[0 0 0 0];
    xd(:,i)=[0 0 0 0];
    xdd(:,i)=inv(m)*((f(i,i)-(c*xd(:,i))-(k*x(:,i))));
   
    end
     
end

    betav= [1/4 1/6];

    for j=1:2
       
    beta = betav(j);
    gamma = 1/2;
   
    end
   
 %%Initial Calculation
a1=(1/(beta*(dt^2)))*m+(gamma/(beta*dt))*c;
a2=(1/(beta*dt))*m +((gamma/beta)-1)*c;
a3=((1/(2*beta))-1)*m +dt*((gamma/(2*beta))-1)*c;
kcap=k+a1;


   for i=1:Z-1
       
    pcap(:,i+1)=f(i+1)+a1*x(:,i)+ a2*xd(:,i)+a3*xdd(:,i);
    x(:,i+1)=(inv(kcap))*pcap(:,i+1);
    xd(:,i+1)=(gamma/(beta*dt))*(x(:,i+1)-x(:,i))+(1-(gamma/beta))*xd(:,i)+dt*(1-(gamma/(2*beta)))*xdd(:,i);
    xdd(:,i+1)=(1/(beta*(dt^2)))*(x(:,i+1)-x(:,i))-(1/(beta*dt))*xd(:,i)-((1/(2*beta))-1)*xdd(:,i);
   
    end

     
    for i=1:n
       
        xddt(i,:)=xdd(i,:)+uddg(i,:);
   
    end

%%Result Plotting
figure(1)
subplot(2,1,1)
plot(t, x(2,:))
grid on
hold on
xlabel('time (sec)')
ylabel('Roof Displacement (m)')
title('Roof Displacement Response (u1)')

subplot(2,1,2)
plot(t,x(1,:))
grid on
hold on
xlabel('time (sec)')
ylabel('First Story Displacement (m)')
title('First Story Displacement Response (u2)')

figure(2)
subplot(2,1,1)
plot(t, xd(2,:))
grid on
hold on
xlabel('time (sec)')
ylabel('Roof Velocity (m/s)')
title('Roof Velocity Response')

subplot(2,1,2)
plot(t, xd(1,:))
grid on
hold on
xlabel('time (sec)')
ylabel('First Story Velocity (m/s)')
title('First Story Velocity Response')

figure(3)
subplot(2,1,1)
plot(t, xddt(2,:))
grid on
hold on
xlabel('time (sec)')
ylabel('Roof total Acceleration (m/s2)')
title('Roof Total Acceleration Response')

subplot(2,1,2)
plot(t, xddt(1,:),'-')
grid on
hold on
xlabel('time (sec)')
ylabel('First Story total Acceleration (m/s2)')
title('First Story Total Acceleration Response')