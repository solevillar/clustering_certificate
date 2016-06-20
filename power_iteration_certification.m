%Author:       Takayuki Iguchi
% Filename:     power_iteration_certification.m
% Last edited:  8 May 2016
% Description:  This function implements the power iteration detector to
%               certify k-means optimality (with k=2) within some small 
%               tolerance epsilon. See [1] for more details. 
% Inputs:       
%               -Phi: 
% 
%               An m x n array of clustered data where m denotes the
%               dimension of the data and n denotes the number of points
%               clustered. The points are sorted according to cluster
%               membership so that 
%                   Phi = [Phi_c1 , Phi_c2];
%               where Phi_cj is an array of all the points assigned to the
%               jth cluster.
% 
%               -nvector: 
% 
%               A 1 x k array denoting the number of points assigned to
%               each cluster. That is, 
%                   nvector = [n_1 , n_2];
%               where n_j denotes the number of points in the jth cluster
%               (corresponding to the clusters in Phi)
% 
%               -epsilon: 
% 
%               A positive number less than 1 (e.g. 10^-6). This parameter 
%               is used to create an upper bound on the probability of 
%               making the error of failing to certify kmeans optimality 
%               when the theoretical result given in the paper should have 
%               certified kmeans optimality. 
% Outputs:
%               -test:
%               
%               Either 0 or 1. A value of 1 denotes certification.
% 
%Documentation:
% 
% [1] Iguchi, Mixon, Peterson, Villar. Probably certifiably correct k-means
%       clustering
% -------------------------------------------------------------------------

function test=power_iteration_certification(Phi,nvector,epsilon)

%initialize stuff
[~,N]   =size(Phi);
na      =nvector(1);
nb      =nvector(2);
Phia    =Phi(:,1:na);
Phib    =Phi(:,na+1:end);

%step 1
nu=zeros(N,1);

for i=1:N
    nu(i)=norm(Phi(:,i))^2;
end
nua     =nu(1:na);
nub     =nu(na+1:end);

%step 4

    %(calculate mu_a)
    temp1   =step2(na,na,ones(na,1),Phia,Phia,nua,nua);
    temp2   =2/na*temp1;
    temp3   =ones(1,na)*temp1;
    temp3   =ones(na,1)*temp3/na^2;
    mua     =(temp3-temp2)/2;

    %(calculate mu_b)
    temp1   =step2(nb,nb,ones(nb,1),Phib,Phib,nub,nub);
    temp2   =2/nb*temp1;
    temp3   =ones(1,nb)*temp1;
    temp3   =ones(nb,1)*temp3/nb^2;
    mub     =(temp3-temp2)/2;

%step 6
Mab1        =step5(ones(nb,1),mua,mub,na,nb,Phia,Phib,nua,nub);
Mba1        =step5(ones(na,1),mub,mua,nb,na,Phib,Phia,nub,nua);
candidate1  =(2*na)/N*min(Mab1);
candidate2  =(2*nb)/N*min(Mba1);
test        =(candidate1<=candidate2);
z           =candidate1*test+candidate2*(1-test);

%step 7
uab         =Mab1-z*N/(2*na)*ones(na,1);
uba         =Mba1-z*N/(2*nb)*ones(nb,1);

%step 8
rhoab       =uab'*ones(na,1);
rhoba       =uba'*ones(nb,1);

%run power method using step 11 to perform matrix vector multiplication
test                =false;
normalized_ones     =ones(1,N)/sqrt(N);
x                   =randn(N,1);
x                   =x/norm(x);
cont                =true;
count               =1;

while cont
    xnew            =step11(x,na,nb,uab,uba,rhoab,rhoba,nu,Phi,N,z);
    ratio           =norm(xnew)/norm(x);
    xnew            =xnew/norm(xnew);
    cos_sq_theta    =(normalized_ones*xnew)^2;
    x               =xnew;
    
    if ratio>z
        cont=false;
    
    elseif cos_sq_theta>1-epsilon
        cont=false;
        test=true;
    
    elseif count>10000*log(1/epsilon)
        fprintf('too many interations\n')
        cont=false;
    
    else
        count=count+1;
    end
end
end
%--------------------------------------------------------------------------
%----------------------------sub funcitons---------------------------------
%--------------------------------------------------------------------------

function output=step2(na,nb,x,Phia,Phib,nua,nub)

%round1
temp1   =ones(1,nb)*x;
temp2   =Phib*x;
temp3   =nub'*x;

%round2
temp1   =nua*temp1;
temp2   =-2*Phia'*temp2;
temp3   =ones(na,1)*temp3;
output  =temp1+temp2+temp3;
end

%--------

function output=step3(nu,Phi,x,N)

%round1
temp1   =nu'*x;
temp2   =Phi*x;
temp3   =ones(1,N)*x;

%round2
temp1   =ones(N,1)*temp1;
temp2   =-2*Phi'*temp2;
temp3   =nu*temp3;
output  =temp1+temp2+temp3;
end

%--------

function output=step5(x,mua,mub,na,nb,Phia,Phib,nua,nub)

%round1
temp1   =mub'*x;
temp2   =ones(1,nb)*x;
temp3   =step2(na,nb,x,Phia,Phib,nua,nub);

%round2
temp1   =ones(na,1)*temp1;
temp2   =mua*temp2;
output  =temp3+temp2+temp1;
end

%-------

function output=step9(x,na,nb,uab,uba,rhoab,rhoba)

xa  =x(1:na);
xb  =x(na+1:end);

%block a
temp1   =uba'*xb;
temp1   =uab*temp1;
temp1   =temp1/rhoba;
    
%block b
temp2   =uab'*xa;
temp2   =uba*temp2;
temp2   =temp2/rhoab;
output  =[temp1;temp2];
end

%-------

function output=step10(x,na,nb)

xa      =x(1:na);
xb      =x(na+1:end);
temp1   =ones(1,na)*xa;
temp1   =ones(na,1)*temp1/na;
temp2   =ones(1,nb)*xb;
temp2   =ones(nb,1)*temp2/nb;
output  =[xa-temp1;xb-temp2];
end

%------

function output=step11(x,na,nb,uab,uba,rhoab,rhoba,nu,Phi,N,z)

temp1   =step10(x,na,nb);
temp2   =step9(temp1,na,nb,uab,uba,rhoab,rhoba);
temp3   =step3(nu,Phi,temp1,N);
temp2   =step10(temp2,na,nb);
temp3   =step10(temp3,na,nb);
temp4   =ones(1,N)*x;
temp4   =z*ones(N,1)*temp4/N;
output  =temp4+(temp3-temp2);
end

