clear all,close all,clc;
%% sostools400
dir1=pwd
cd('C:\Users\Artun\Documents\MATLAB\SOSTOOLS400')
run('addsostools.m')
cd(dir1)
clear all,close all,clc;
%% this is for "SDPT3"
current_dir=pwd;
cd('C:\Users\Artun\Documents\MATLAB\sdpt3');
run('install_sdpt3.m')
cd(current_dir);
clear all,close all,clc;
%% mosek 2015a
addpath(genpath('C:\Program Files\Mosek\9.3\toolbox\R2015a'))
%% yalmip
addpath(genpath('C:\Users\Artun\Documents\MATLAB\YALMIP-master'))
% "clear classes"
% run('yalmiptest.m')

%% BEFORE RUNNING THE CODE ADD "YALMIP" AND "SDPT3" libraries
%     set(findall(gcf,'type','line'),'linewidth',[1]);
% yalmip('clear');
%%

%% add MOSEK
%% add YALMIP

clear all,close all,clc;

%% step-1 compute A(t)
syms x1 x2 real

f1=x2;
f2=-x1-2*x2+0.25*(x1^2)*x2;

A_t=jacobian([f1;f2],[x1;x2]);
syms x1hat x2hat real
A_t=subs(A_t,[x1 x2],[x1hat x2hat])
%% step-2 decide Q,P(0),R
Q=eye(2)    % nxn
R=eye(1)    % pxp
C=[1,0]
B=[0;1]
%% step-3 Pdot-dyns
syms p11 p12 p22 real
P=[p11,p12;p12,p22]
Pdot=[A_t*P]+[A_t*P]'+Q-P*C'*inv(R)*C*P

Pdot(1,1)   % p11_dot
Pdot(1,2)   % p12_dot
Pdot(2,2)   % p22_dot
%% step-4 compute H(t)
H_t=P*C'*inv(R)
%% step-5 explicity write down xhat_dot-dyns
H_t*[x1-x1hat]

[subs(f1,[x1 x2],[x1hat x2hat]);subs(f2,[x1 x2],[x1hat x2hat])]+H_t*[x1-x1hat]


%%

fig1=figure(1);fig1.Color=[1,1,1];
ax1=axes('Parent',fig1);
    set(0,'CurrentFigure',fig1);
    set(fig1,'currentaxes',ax1);
for ii=1:1:1
tspan=[0:0.01:10]; x0=randn(7,1)*1; x0(5)=1;x0(6)=0;x0(7)=1;
wt=tspan;
% f=randi([1,20],1,1); w=sin(2*pi*f*tspan);
w=square(tspan);
[t,x]=ode45(@(t,x) odefcn(t,x,wt,w),tspan,x0);
x1_vec=x(:,1);x2_vec=x(:,2);
x1hat_vec=x(:,3);x2hat_vec=x(:,4);
plot(t,x1_vec,'r-','LineWidth',[2],"Parent",ax1);hold on;
plot(t,x1hat_vec,'g-','LineWidth',[2],"Parent",ax1);hold on;
hold on;yline(1);yline(-1);
end
function xdot=odefcn(t,x,wt,w)
w=interp1(wt,w,t);
xdot=zeros(7,1);
x1=x(1);x2=x(2);
x1hat=x(3);x2hat=x(4);
p11=x(5);p12=x(6);p22=x(7);

f1=x2;
f2=-x1-2*x2+0.25*(x1)^2*(x2)+0.2*sin(2*t);
x1dot=f1;
x2dot=f2;
y=x1;
x1hat_dot=x2hat + p11*(x1 - x1hat);
x2hat_dot=- 2*x2hat - x1hat + (x1hat^2*x2hat)/4+0.2*sin(2*t)+p12*(x1 - x1hat);
%
p11dot=- p11^2 + 2*p12 + 1;
p12dot=p22 - p11*p12 + p12*(x1hat^2/4 - 2) + p11*((x1hat*x2hat)/2 - 1);
p22dot=2*p22*(x1hat^2/4 - 2) + 2*p12*((x1hat*x2hat)/2 - 1) - p12^2 + 1;


xdot(1)=f1;
xdot(2)=f2;
xdot(3)=x1hat_dot;
xdot(4)=x2hat_dot;
xdot(5)=p11dot;
xdot(6)=p12dot;
xdot(7)=p22dot;
end






















%