clc
clear all
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
%CP1(:,31:35)=[];
%%
CP1(:,14)=[];
CP1(:,32:36)=[];

%%
%Canot bring France If I am to include oil and commodity , France causes
%singularity in both.

%%
%This version executes. We can keep UK, but we need to keep either oil or
%commodity. Meaning UK is causing singularity to commodity.
% CP1(:,32:35)=[];
% CP1(:,32)=[];
%%
%This version executes. Meaning, We can keep UK but we need to remove
%commodity. Which means UK causing singularity to oil.
% CP1(:,32:35)=[];
% CP1(:,33)=[];
%% Do not forget to include UK. While doing separate testing include UK. Thats important.
%%
%Remember when u omit number 1 column, number 2 becomes number 1. 
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
windowSize=100;
%%
% creating empty array to store results
contToOthers=[];
contFromOthers=[];
%%
%for J = 1:windowSize: (nRow-windowSize) % size by size

    y = filt;     
    z = bsxfun(@rdivide,bsxfun(@minus,y,mean(y)),std(y));

    %disp('    Mean      Median    Max       Min       Std.dev.  Skew      Kurt')
    %disp( [ mean(y)' median(y)' max(y)' min(y)' std(y)' mean(z.^3)' mean(z.^4)'] )
    x = [ones(length(y)-2,1)   trimr(y,1,1)   trimr(y,0,2) ];
    b = x\trimr(y,2,0);
 
    %% Taking total number of rows and columns 
    [rcnt, ccnt] = size(b);
%%
    mu  = trimr(b,0,rcnt-1)';  % Vector of intercepts                                              
    phi1 = trimr(b,1,ccnt)'; % Lag 1 parameter estimates   
    phi2 = trimr(b,ccnt+1,0)'; % Lag 2 parameter estimates     

    %% Generate VMA (non-orthogonalized) for horizons 1 to 10 
    si1  = eye(size(y,2));
    si2  = phi1;
    si3  = phi1*si2 + phi2;
    si4  = phi1*si3 + phi2*si2;
    si5  = phi1*si4 + phi2*si3;
    si6  = phi1*si5 + phi2*si4;
    si7  = phi1*si6 + phi2*si5;
    si8  = phi1*si7 + phi2*si6;
    si9  = phi1*si8 + phi2*si7;
    si10  = phi1*si9 + phi2*si8;
%% orthogonalization
    v  = trimr(y,2,0) - x*b;   % VAR residuals      
    vc = v'*v/length(v);
    d  = diag(vc);
%    s  = chol(vc);% this is what makes the shocks structural or linear or orthogonalized
%% non-orthogonalized GIRF (meaning they are not effected by variable ordering)
    ir1  = si1*vc /(sqrt(d).*eye(size(vc,2)));
    ir2  = si2*vc /(sqrt(d).*eye(size(vc,2)));
    ir3  = si3*vc /(sqrt(d).*eye(size(vc,2)));
    ir4  = si4*vc /(sqrt(d).*eye(size(vc,2)));
    ir5  = si5*vc /(sqrt(d).*eye(size(vc,2)));
    ir6  = si6*vc /(sqrt(d).*eye(size(vc,2)));
    ir7  = si7*vc /(sqrt(d).*eye(size(vc,2)));
    ir8  = si8*vc /(sqrt(d).*eye(size(vc,2)));
    ir9  = si9*vc /(sqrt(d).*eye(size(vc,2)));
    ir10 = si10*vc /(sqrt(d).*eye(size(vc,2)));
%%
impulse = ir10';
impulse=impulse(:);
[m,n]=size(v);
impulse1=repmat(impulse,13);
u=v(:);
impulse1=impulse1(:);
fixing= length(impulse1)-length(u);
impulse1=trimr(impulse1,0,fixing);
impulse1= impulse1(:);
u=reshapeg(u,[],31);% Change 31 as per number of variables
impulse=reshapeg(impulse1,[],31); % change 31 as per number of variabls.
GHDTO = u* impulse'; 
GHDto= bsxfun(@rdivide,GHDTO,sum(GHDTO));
scalecheckto= sum(GHDto);% if columns are summing upto 1 meaning, shock j gives to i
GHDFROM= u'.*impulse';

GHDfrom= bsxfun(@rdivide,GHDFROM,sum(GHDFROM,2));
scalecheckfrom=sum(GHDfrom,2);% If we get 1's in columns, meaning tnis is shock i receives from j

%% Adding somethimg new to the code
TO= GHDto(:,1:31);
FROM= GHDfrom(:,1:31);


%% without scaling to 1 perfectly works, use this one
plot(movmean(GHDTO(:,1),250))
hold on
plot(movmean(GHDFROM(1,:),250))
hold off
%% with scaling to one


% %%
% plot(movmean(GHD(:,1),250))
% hold on
% plot(movmean(GHD(1,:),250))
% hold off
% 
% GHDFROM= u'.*impulse';
% % %GHD is shocks to others
% GHDTO=bsxfun(@rdivide,GHDTO,sum(GHDTO));
% % save GHDTO.mat;
% % %GHD is shocks received from others
% GHDFROM=bsxfun(@rdivide,GHDFROM,sum(GHDFROM));
% % save GHDFROM.mat;
% % dlmwrite('GHD_shock_to_others.csv',GHDTO);
% % dlmwrite('GHD_shock_from_others.csv',GHDFROM);
% 
% plot(movmean(GHDTO(:,1),250))
% hold on
% plot(movmean(GHDFROM(1,:),250))
% hold off
% 
% 
  