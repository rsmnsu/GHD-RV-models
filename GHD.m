clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
CP1(:,32:36)=[];

% ind=GHDTO(:,3);
% mys=GHDTO(:,5);
% sgp=GHDTO(:,7);
% kor=GHDTO(:,9);
% phl=GHDTO(:,8);
% tha=GHDTO(:,11);4

CP1=CP1(:,[3 5 7 9 8 11 28 29 31 17 30 2 20 21 22 23 24 25 1 19 26 27 4 6 16 18 15 10 12 13 14]);
%%
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
y = filt; 
 p = 2;      % Order of VAR lags
    q = 2;     % Order of VMA lags  
    %%
    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1) ;
    nconst = size(ylag,2);
    for K = 1:p
        ylag = [ trimr(ylag,1,0) trimr(y,0,K) ];
    end
    %%
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;   
    omegav = v'*v/size(v,1);
%%
ncol = size(v,2);

u(:,1) = v(:,1);
a = zeros(ncol, ncol);
for K = 2 : size(v,2)
  a(1:K-1,K) = v(:,1:K-1)\v(:,K);
  u(:,K) = v(:,K) - v(:,1:K-1)*a(1:K-1,K);
end

b = diag(ones(1,ncol)) - tril(a.');

residual=u;
%%
% Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(u,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
%     s=chol(omegav);
  %%
    % Construct C(L) matrices (vector moving average) 
       c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ])

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end

       tmp = ss';
       c = [ c  tmp(:) ];

    end

%%    
   ar0=a1;
   ar0={ar0};
   ma0=tmp;
   ma0={ma0};
   InnovCov=omegav; 
   GIRF = armairf(ar0,ma0,'Method','generalized','InnovCov',InnovCov,'NumObs',10);
   
    %% The old stuff begins here, lets see
 
impulse=GIRF;
impulse1=permute(impulse, [1 3 2]);
impulse2=reshape(impulse1, [], size(impulse, 2),1);
impulse2=impulse(:);
impulse3= repmat(impulse2,10);

fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
u=reshapeg(u,[],31);
impulse5=reshapeg(impulse4,[],31);
%% experiment with To and From
GHDtransmission = u* impulse5'; %plot(movmean(HD(:,1),250)) works
GHDTO= GHDtransmission(:,1:31);
GHDvulnerability= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
GHDFROM=GHDvulnerability';

GHDCTO=GHDTO(:,1:24);
GHDCFO= GHDFROM(:,1:24);
save GHDCTO.mat GHDCTO;
save GHDCFO.mat GHDCFO;
GHDtransmission= movmean(GHDTO, 250);
GHDvulnerability= movmean(GHDFROM,250);


%% Plots

% %% AC
% %India
% subplot (2,3,1)
% plot(GHDtransmission(:,1))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,1))
% hold off
% %Malaysia
% subplot (2,3,2)
% plot(GHDtransmission(:,2))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,2))
% hold off
% %Singapore
% subplot (2,3,3)
% plot(GHDtransmission(:,3))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,3))
% hold off
% %S.korea
% subplot(2,3,4)
% plot(GHDtransmission(:,4))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,4))
% hold off
% %Phillipines
% subplot(2,3,5)
% plot(GHDtransmission(:,5))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,5))
% hold off
% %Thailand
% subplot(2,3,6)
% plot(GHDtransmission(:,6))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,6))
% hold off

%% Plot with functions
%% AC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(GHDtransmission(:,1),'Parent',subplot1,'DisplayName','Transmission');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(GHDvulnerability(:,1),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('India');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(GHDtransmission(:,2),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(GHDvulnerability(:,2),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Malaysia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(GHDtransmission(:,3),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(GHDvulnerability(:,3),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Singapore');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(GHDtransmission(:,4),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(GHDvulnerability(:,4),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('South Korea');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(GHDtransmission(:,5),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(GHDvulnerability(:,5),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('The Phillipines');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(GHDtransmission(:,6),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(GHDvulnerability(:,6),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Thailand');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');


% %% EC
% % Germany
% subplot (2,3,1)
% plot(GHDtransmission(:,7))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,7))
% hold off
% % Chile
% subplot (2,3,2)
% plot(GHDtransmission(:,8))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,8))
% hold off
% %France
% subplot (2,3,3)
% plot(GHDtransmission(:,9))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,9))
% hold off
% % China
% subplot(2,3,4)
% plot(GHDtransmission(:,10))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,10))
% hold off
% % UK
% subplot(2,3,5)
% plot(GHDtransmission(:,11))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,11))
% hold off
% % Australia
% subplot(2,3,6)
% plot(GHDtransmission(:,12))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,12))
% hold off

%% EC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(GHDtransmission(:,7),'Parent',subplot1,'DisplayName','Transmission');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(GHDvulnerability(:,7),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Germany');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(GHDtransmission(:,8),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(GHDvulnerability(:,8),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Chile');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(GHDtransmission(:,9),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(GHDvulnerability(:,9),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('France');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(GHDtransmission(:,10),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(GHDvulnerability(:,10),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('China');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(GHDtransmission(:,11),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(GHDvulnerability(:,11),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('UK');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(GHDtransmission(:,12),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(GHDvulnerability(:,12),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Australia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
%% GC
% % Greece
% subplot (2,3,1)
% plot(GHDtransmission(:,13))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,13))
% hold off
% % Portugal
% subplot (2,3,2)
% plot(GHDtransmission(:,14))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,14))
% hold off
% % Ireland
% subplot (2,3,3)
% plot(GHDtransmission(:,15))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,15))
% hold off
% % Belgium
% subplot(2,3,4)
% plot(GHDtransmission(:,16))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,16))
% hold off
% % Croatia
% subplot(2,3,5)
% plot(GHDtransmission(:,17))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,17))
% hold off
% % Austria
% subplot(2,3,6)
% plot(GHDtransmission(:,18))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,18))
% hold off


%% GC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(GHDtransmission(:,13),'Parent',subplot1,'DisplayName','Transmission');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(GHDvulnerability(:,13),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Greece');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(GHDtransmission(:,14),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(GHDvulnerability(:,14),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Portugal');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(GHDtransmission(:,15),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(GHDvulnerability(:,15),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Ireland');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(GHDtransmission(:,16),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(GHDvulnerability(:,16),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Belgium');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(GHDtransmission(:,17),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(GHDvulnerability(:,17),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Croatia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(GHDtransmission(:,18),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(GHDvulnerability(:,18),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Austria');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');

% %% OED
% % The USA
% subplot (2,3,1)
% plot(GHDtransmission(:,19))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,19))
% hold off
% % Canada
% subplot (2,3,2)
% plot(GHDtransmission(:,20))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,20))
% hold off
% % Russia
% subplot (2,3,3)
% plot(GHDtransmission(:,21))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,21))
% hold off
% % Norway
% subplot(2,3,4)
% plot(GHDtransmission(:,22))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,22))
% hold off
% % Japan
% subplot(2,3,5)
% plot(GHDtransmission(:,23))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,23))
% hold off
% % New Zealand
% subplot(2,3,6)
% plot(GHDtransmission(:,24))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,24))
% hold off


%% OED
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(GHDtransmission(:,19),'Parent',subplot1,'DisplayName','Transmission');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(GHDvulnerability(:,19),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('The USA');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(GHDtransmission(:,20),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(GHDvulnerability(:,20),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Canada');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(GHDtransmission(:,21),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(GHDvulnerability(:,21),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Russia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(GHDtransmission(:,22),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(GHDvulnerability(:,22),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Norway');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(GHDtransmission(:,23),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(GHDvulnerability(:,23),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Japan');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(GHDtransmission(:,24),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(GHDvulnerability(:,24),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('New Zealand');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');

% %% OEE
% % Saudi Arabia
% subplot (2,3,1)
% plot(GHDtransmission(:,25))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,25))
% hold off
% % Israel
% subplot (2,3,2)
% plot(GHDtransmission(:,26))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,26))
% hold off
% % Iraq
% subplot (2,3,3)
% plot(GHDtransmission(:,27))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,27))
% hold off
% % Sri Lanka
% subplot(2,3,4)
% plot(GHDtransmission(:,28))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,28))
% hold off
% % Nigeria
% subplot(2,3,5)
% plot(GHDtransmission(:,29))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,29))
% hold off
% % Venezuela
% subplot(2,3,6)
% plot(GHDtransmission(:,30))
% hold on
% yyaxis right
% plot(GHDvulnerability(:,30))
% hold off
%% OEE
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(GHDtransmission(:,25),'Parent',subplot1,'DisplayName','Transmission');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(GHDvulnerability(:,25),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Saudi Arabia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(GHDtransmission(:,26),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(GHDvulnerability(:,26),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Israel');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(GHDtransmission(:,27),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(GHDvulnerability(:,27),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Iraq');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(GHDtransmission(:,31),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(GHDvulnerability(:,31),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Kuwait');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(GHDtransmission(:,29),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(GHDvulnerability(:,29),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Nigeria');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(GHDtransmission(:,30),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(GHDvulnerability(:,30),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Venezuela');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
%% Save 
% HD=zscore(HD);
% HD1=zscore(HD1);
save GHDtransmission.mat GHDTO ;
save GHDvulnerability.mat GHDFROM ;



%% Maps: begins here is ANN

X=GHDtransmission;
N=40;
R=size(X,1);
pos=floor(linspace(1,R,N+1));
pos(end)=R+1;
sizes=diff(pos);
subsets=mat2cell(X,sizes,size(X,2));
A1=subsets{1,1};
A2=subsets{2,1};
A3=subsets{3,1};
A4=subsets{4,1};
A5=subsets{5,1};
A6=subsets{6,1};
A7=subsets{7,1};
A8=subsets{8,1};
A9=subsets{9,1};
A10=subsets{10,1};
A11=subsets{11,1};
A12=subsets{12,1};
A13=subsets{13,1};
A14=subsets{14,1};
A15=subsets{15,1};
A16=subsets{16,1};
A17=subsets{17,1};
A18=subsets{18,1};
A19=subsets{19,1};
A20=subsets{20,1};
A21=subsets{21,1};
A22=subsets{22,1};
A23=subsets{23,1};
A24=subsets{24,1};
A25=subsets{25,1};
A26=subsets{26,1};
A27=subsets{27,1};
A28=subsets{28,1};
A29=subsets{29,1};
A30=subsets{30,1};
A31=subsets{31,1};
A32=subsets{32,1};
A33=subsets{33,1};
A34=subsets{34,1};
A35=subsets{35,1};
A36=subsets{36,1};
A37=subsets{37,1};
A38=subsets{38,1};
A39=subsets{39,1};
A40=subsets{40,1};
%% Full network
inputs = GHDtransmission;

% Create a Self-Organizing Map
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);

% Train the Network
[net,tr] = train(net,inputs);

% Test the Network
outputs = net(inputs);

% View the Network
view(net)

% Plots
% Uncomment these lines to enable various plots.
% figure, plotsomtop(net)
% figure, plotsomnc(net)
fig= plotsomnd(net);
savefig('GHDtransmissionfull.fig');
saveas (fig,'GHDtransmissionfull.tif');
% figure, plotsomplanes(net)
% figure, plotsomhits(net,inputs)
% figure, plotsompos(net,inputs)

%% A1
inputs = A1;

% Create a Self-Organizing Map
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);

% Train the Network
[net,tr] = train(net,inputs);

% Test the Network
outputs = net(inputs);

% View the Network
view(net)

% Plots
% Uncomment these lines to enable various plots.
% figure, plotsomtop(net)
% figure, plotsomnc(net)
fig= plotsomnd(net);
savefig('A1.fig');
saveas (fig,'A1.tif');
% figure, plotsomplanes(net)
% figure, plotsomhits(net,inputs)
% figure, plotsompos(net,inputs)

%% A2
inputs = A2;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A2.fig');
saveas (fig,'A2.tif');
%% A3
inputs = A3;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A3.fig');
saveas (fig,'A3.tif');
%% A4
inputs = A4;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A4.fig');
saveas (fig,'A4.tif');
%% A5
inputs = A5;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A5.fig');
saveas (fig,'A5.tif');
%% A6
inputs = A6;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A6.fig');
saveas (fig,'A6.tif');
%% A7
inputs = A7;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A7.fig');
saveas (fig,'A7.tif');
%% A8
inputs = A8;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A8.fig');
saveas (fig,'A8.tif');
%% A9
inputs = A9;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A9.fig');
saveas (fig,'A9.tif');
%% A10
inputs = A10;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A10.fig');
saveas (fig,'A10.tif');
%% A11
inputs = A11;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A11.fig');
saveas (fig,'A11.tif');
%% A12
inputs = A12;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A12.fig');
saveas (fig,'A12.tif');
%% A13
inputs = A13;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A13.fig');
saveas (fig,'A13.tif');
%% A14
inputs = A14;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A14.fig');
saveas (fig,'A14.tif');
%% A15
inputs = A15;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A15.fig');
saveas (fig,'A15.tif');
%% A16
inputs = A16;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A16.fig');
saveas (fig,'A16.tif');
%% A17
inputs = A17;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A17.fig');
saveas (fig,'A17.tif');
%% A18
inputs = A18;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A18.fig');
saveas (fig,'A18.tif');
%% A19
inputs = A19;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A19.fig');
saveas (fig,'A19.tif');
%% A20
inputs = A20;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A20.fig');
saveas (fig,'A20.tif');
%% A21
inputs = A21;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A21.fig');
saveas (fig,'A21.tif');
%% A22
inputs = A22;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A22.fig');
saveas (fig,'A22.tif');
%% A23
inputs = A23;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A23.fig');
saveas (fig,'A23.tif');
%% A24
inputs = A24;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A24.fig');
saveas (fig,'A24.tif');
%% A25
inputs = A25;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A25.fig');
saveas (fig,'A25.tif');
%% A26
inputs = A26;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A26.fig');
saveas (fig,'A26.tif');
%% A27
inputs = A27;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A27.fig');
saveas (fig,'A27.tif');
%% A28
inputs = A28;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A28.fig');
saveas (fig,'A28.tif');
%% A29
inputs = A29;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A29.fig');
saveas(fig,'A29.tif');
%% A30
inputs = A30;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A30.fig');
saveas(fig,'A30.tif');
%% A31
inputs = A31;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A31.fig');
saveas(fig,'A31.tif');
%% A32
inputs = A32;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A32.fig');
saveas(fig,'A32.tif');
%% A33
inputs = A33;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A33.fig');
saveas(fig,'A33.tif');
%% A34
inputs = A34;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A34.fig');
saveas(fig,'A34.tif');
%% A35
inputs = A35;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A35.fig');
saveas(fig,'A35.tif');
%% A36
inputs = A36;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A36.fig');
saveas(fig,'A36.tif');
%% A37
inputs = A37;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A37.fig');
saveas(fig,'A37.tif');
%% A38
inputs = A38;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A38.fig');
saveas(fig,'A38.tif');
%% A39
inputs = A39;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A39.fig');
saveas(fig,'A39.tif');
%% A40
inputs = A40;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A40.fig');
saveas(fig,'A40.tif');



