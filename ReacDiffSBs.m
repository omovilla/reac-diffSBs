clc
clear all

%Parameters
dx=0.05;
xmin=-1;
xmax=1;
xmesh = [xmin:dx:xmax];
dt=0.2;
T=10;
tspan = [0:dt:T];
tspanb = [T:-dt:0];
[X,Y]=meshgrid(tspan,xmesh);
sigma=0.1; 
N=length(xmesh)+1;
iter=5; % number of iterations


% defining a Gaussian distribution
U0=@(t,x) 1/2*(x).^2/0.1;
Z=sum(exp(-U0(0,xmesh))).*dx;
q0=@(x)exp(-U0(0,x))/Z;


% defining initial and final conditions (Schrodinger bridge conditions)
p0=@(x)q0(x); % initial condition on the primary space
pf=@(x)0.5.*q0(x); % final condition on the primary space
M=1000000; % total mass of the system
m0=dx.*sum(p0(xmesh)); % mass on the primary space at the beginning
mf=dx.*sum(pf(xmesh)); % mass on the primary space at the end
mu0=M-m0; % initial condition on reservoir
muf=M-mf; % final condition on reservoir


% Initializing iterations
phi=ones(length(tspan),length(xmesh)); %values of phi(t,x) for the first iteration
psi=ones(length(tspan),1); %values of psi(t,x) for the first iteration

hatPhi0=[p0(xmesh),mu0]; %initial condition for the forward equation

for i=1:iter

% Solving the forward equation
[t,hatPhi]=ode45(@(t,y)odehat(t,y,sigma,M,N,dx,tspan,phi,psi),tspan,hatPhi0);

% Finding the final condition for the backward equation
PhiT=[pf(xmesh),muf]./hatPhi(end,:);

% Solving the backward equation
[t,Phi]=ode45(@(t,y)odeback(t,y,sigma,M,N,dx,tspan,hatPhi),tspanb,PhiT);

% Writting the solution to the backward eq forward in time
for j=1:length(t)
    phi(j,:)=Phi(end-j+1,[1:N-1]);
    psi(j,:)=Phi(end-j+1,end);
end

% Finding the new initial condition for the forward equation
hatPhi0=[p0(xmesh),mu0]./Phi(end,:);

i
end

%Plot of primary space density
figure
surf(xmesh,tspan,phi.*hatPhi(:,[1:N-1]),'EdgeColor', 'none')
alpha 0.9
colormap(pink)
hold on
plot3(xmesh,tspan(1)*ones(1,length(xmesh)),p0(xmesh),'LineWidth',3)
plot3(xmesh,tspan(end)*ones(1,length(xmesh)),pf(xmesh),'LineWidth',3)

%Plot of reservoir density
figure
plot(tspan,psi.*hatPhi(:,end),'LineWidth',3)

%%%%%%


%Solves the FORWARD PDE through the method of lines
function dydt = odehat(t,y,sigma,M,N,dx,tspan,phi,psi)

%each ODE corresponds to a different point in the primary space (the last one corresponds to the reservor)
dydt = zeros(N,1); 

%Initializing the spatial integrals needed for the reservoir ODE
intA=0;
intD=0;

%Loop to define the ODES on the primary space
for i=1:N-1

    %Laplacian with no flux boundary    
    if i==1
       dxxphihat=(y(i+1)-2.*y(i)+y(i))./dx^2;
    elseif i==N-1;
       dxxphihat=(y(i)-2.*y(i)+y(i-1))./dx^2;
    else
       dxxphihat=(y(i+1)-2.*y(i)+y(i-1))./dx^2;
    end

    %phi and psi need to be functions for the ODE solver
    phifun=@(t)interp1(tspan,phi(:,i),t);
    psifun=@(t)interp1(tspan,psi,t);
    
    %one-time marginal at each gridpoint in space as a function of time
    p=phifun(t).*y(i);
    
    %PRIOR PARAMETERS (since they depend on p, these need to be changed manualy in the forward and backward equations)
    A=p./M;% Parameters corresponding to the Fisher equation with r=1
    D=p;%
    dA=1./M;% Derivative of A with respect to p
    dD=1;% Derivative of D with respect to p

    %ODEs of the primary space
    dydt(i) = +1/2.*sigma.^2.*dxxphihat-y(i).*(D+(1-psifun(t)./phifun(t)).*p.*dD)+y(N).*(A+(1-psifun(t)./phifun(t)).*p.*dA);
    
    %Integrals for the ODE in the reservoir
    intA=A.*dx+intA;
    intD=y(i).*D.*dx+intD;
end

%reservoir ODE
dydt(N) = -y(N).*intA +intD;
end


%Solves the BACKWARD PDE through the method of lines
function dydt = odeback(t,y,sigma,M,N,dx,tspan,hatPhi)

%each ODE corresponds to a different point in the primary space (the last one corresponds to the reservor)
dydt = zeros(N,1); 

%Initializing the spatial integrals needed for the reservoir ODE
intA=0;

%Loop to define the ODES on the primary space
for i=1:N-1

    %Laplacian with no flux boundary    
    if i==1
      dxxphi=(y(i+1)-2.*y(i)+y(i))./dx^2;
    elseif i==N-1;
      dxxphi=(y(i)-2.*y(i)+y(i-1))./dx^2;
    else
      dxxphi=(y(i+1)-2.*y(i)+y(i-1))./dx^2;
    end

    %phi and psi need to be functions for the ODE solver
    hatpsifun=@(t)interp1(tspan,hatPhi(:,end),t);   
    hatphifun=@(t)interp1(tspan,hatPhi(:,i),t);

    %one-time marginal at each gridpoint in space as a function of time
    p=hatphifun(t).*y(i);

    %PRIOR PARAMETERS (since they depend on p, these need to be changed manualy in the forward and backward equations)
    A=p./M;% Parameters corresponding to the Fisher equation with r=1
    D=p;%
    dA=1./M;% Derivative of A with respect to p
    dD=1;% Derivative of D with respect to p

    %ODEs of the primary space
    dydt(i) = -1/2.*sigma.^2.*dxxphi+(y(i)-y(N)).*(D+p.*dD)+hatpsifun(t)./hatphifun(t).*(y(N)-y(i)).*p.*dA;

    %Integrals for the ODE in the reservoir
    intA=(y(N)-y(i)).*A.*dx+intA;
end

%reservoir ODE
dydt(N) = intA;
end



