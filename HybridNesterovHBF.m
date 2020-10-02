%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HybridNessterovHBF.m
%--------------------------------------------------------------------------
% Project: Uniting Nesterov's accelerated gradient descent globally with
% heavy ball locally. Nonstrongly convex version.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2020 1:07:00

clear all

% global variables
global delta M gamma lambda c_0 c_10 r tauMin tauMax tauMed c 

%%%%%%%%%%%%%%%%%%%%%
% setting the globals
%%%%%%%%%%%%%%%%%%%%%
setMinima();

% Nesterov constants
M = 2;

% Heavy Ball constants
gamma = 2/3; % 
lambda = 40; % 

c_0 = 400; % \mathcal{U}_0 
c_10 = 111; % \mathcal{T}_{1,0} 

Root = (-3 + sqrt(29))/2;

tauMin = 3/(3/((Root+2) - (Root - 1)/(Root+1)));

c = 0.25;
deltaMed = 4236; 
r = 51; 

delta = 0.5;

%%%%%%%%%%%%%%%%%%%%%
% setting the locals
%%%%%%%%%%%%%%%%%%%%%

deltaVecUniting = [0,0,0];
deltaVec = [0,0,0];

% initial conditions
z1_0 = 50; 
z2_0 = 0;
z2_00 = 50;
q_0 = 0;
tau_0 = 0;
tauPN_0 = tauMin; 

tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin; 
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN=[0 400];
JSPAN = [0 20000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate
[tNest,jNest,xNest] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options);

deltaVecNest = timeToConv(xNest,tNest);

[tHBF,jHBF,xHBF] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN,JSPAN,rule,options);

deltaVecHBF = timeToConv(xHBF,tHBF);

[tUniting,jUniting,xUniting] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options);

deltaVecUniting = timeToConv(xUniting,tUniting);

[t,j,x] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options);

deltaVec = timeToConv(x,t);

minarc = min([length(x),length(xUniting)]);
ta = [tUniting(1:minarc),t(1:minarc)];
ja = [jUniting(1:minarc),j(1:minarc)];
xa = [xUniting(1:minarc,1),x(1:minarc,1)];
xb = [xUniting(1:minarc,2),x(1:minarc,2)];

figure(1) 
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 3;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 3;
subplot(2,1,1), plotHarc(ta,ja,xa,[],modificatorF,modificatorJ);
hold on
plot(deltaVec(3),deltaVec(1),'k.','MarkerSize', 20)
strDelta = [num2str(deltaVec(3)), 's'];
text(deltaVec(3),deltaVec(1),strDelta,'HorizontalAlignment','left','VerticalAlignment','top');
plot(deltaVecUniting(3),deltaVecUniting(1),'k.','MarkerSize', 20)
strDeltaUniting = [num2str(deltaVecUniting(3)), 's'];
text(deltaVecUniting(3),deltaVecUniting(1),strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([0 10 -50 70])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
hold off
subplot(2,1,2), plotHarc(ta,ja,xb,[],modificatorF,modificatorJ);
hold on
plot(deltaVec(3),deltaVec(2),'k.','MarkerSize', 20)
plot(deltaVecUniting(3),deltaVecUniting(2),'k.','MarkerSize', 20)
axis([0 10 -100 60])
grid on
ylabel('z_2','Fontsize',16)
xlabel('t','Fontsize',16)
hold off
saveas(gcf,'Plots\ComparisonPlotsNSC','png')

minarc = min([length(xUniting),length(x),length(xHBF),length(xNest)]);
tc = [tUniting(1:minarc),t(1:minarc),tHBF(1:minarc),tNest(1:minarc)];
jc = [jUniting(1:minarc),j(1:minarc),jHBF(1:minarc),jNest(1:minarc)];
xc = [xUniting(1:minarc,1),x(1:minarc,1),xHBF(1:minarc,1),xNest(1:minarc,1)];
xd = [xUniting(1:minarc,2),x(1:minarc,2),xHBF(1:minarc,2),xNest(1:minarc,2)];

figure(2)
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 3;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 3;
subplot(2,1,1), plotHarc(tc,jc,xc,[],modificatorF,modificatorJ);
hold on
plot(deltaVec(3),deltaVec(1),'k.','MarkerSize', 20)
strDelta = [num2str(deltaVec(3)), 's'];
text(deltaVec(3),deltaVec(1),strDelta,'HorizontalAlignment','left','VerticalAlignment','bottom');
plot(deltaVecUniting(3),deltaVecUniting(1),'k.','MarkerSize', 20)
strDeltaUniting = [num2str(deltaVecUniting(3)), 's'];
text(deltaVecUniting(3),deltaVecUniting(1),strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom');
plot(deltaVecNest(3),deltaVecNest(1),'k.','MarkerSize', 20)
strDeltaNest = [num2str(deltaVecNest(3)), 's'];
text(deltaVecNest(3),deltaVecNest(1),strDeltaNest,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([0 15 -20 70])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
axes('Position',[0.7 0.78 0.15 0.08])
box on
hold on
plot(tHBF,xHBF(:,1),'LineWidth',3)
plot(deltaVecHBF(3),deltaVecHBF(1),'k.','MarkerSize', 20)
strDeltaHBF = [num2str(deltaVecHBF(3)), 's'];
text(deltaVecHBF(3),deltaVecHBF(1),strDeltaHBF,'HorizontalAlignment','left','VerticalAlignment','bottom');
hold off
set(gca,'xtick',[])
set(gca,'ytick',[])
axis([0 200 -20 70])
hold off
subplot(2,1,2), plotHarc(tc,jc,xd,[],modificatorF,modificatorJ);
hold on
plot(deltaVec(3),deltaVec(2),'k.','MarkerSize', 20)
plot(deltaVecHBF(3),deltaVecHBF(2),'k.','MarkerSize', 20)
plot(deltaVecNest(3),deltaVecNest(2),'k.','MarkerSize', 20)
plot(deltaVecUniting(3),deltaVecUniting(2),'k.','MarkerSize', 20)
axis([0 15 -50 70])
grid on
ylabel('z_2','Fontsize',16)
xlabel('t','Fontsize',16)
hold off
saveas(gcf,'Plots\ComparisonPlots2','png')
saveas(gcf,'Plots\ComparisonPlots2','epsc')
