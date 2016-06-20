% Author:       Takayuki Iguchi
% Filename:     kmeans_certificate_balls.m
% Last edited:  14 May 2016
% Description:  This function generates points uniformly from k balls in an
%               m dimensional space at a desired ball center distance. 
%               The function then implements much of the procedure outlined
%               in Section 3.2 of [1] to test the sufficient condition for
%               k-means optimality in Theorem 7 in [1]. See [1] for more 
%               details. 
% Inputs:       
%               -Delta:
%               
%               The inter-cluster center distance. 
% 
%               -nvector:
% 
%               A 1 x k array where the ith entry corresponds to the number
%               of points in the ith ball.
% 
%               m:
% 
%               The dimension of the points to be generated.
% Outputs:
%               -Theorem7check:
% 
%               Either a 0 or 1 where a value of 1 indicates certification 
%               k-means optimality
%Documentation:
% 
% [1] Iguchi, Mixon, Peterson, Villar. Probably certifiably correct k-means
%       clustering
% -------------------------------------------------------------------------

function [Theorem7check]=kmeans_certificate_balls(Delta,nvector,m)

%initialize--------------
k       =max(size(nvector));
n_sum   =nvector(1)*ones(1,k);

for i=2:k
    n_sum(i)=nvector(i)+n_sum(i-1);
end
n_sum   =horzcat(0,n_sum);
c       =eye(k,m)*Delta/sqrt(2);

% make_balls4calc--------------
s       =size(c);
M       =s(2);%dimension
c       =c';
s       =size(nvector);
k       =max(s);%number of clusters
N       =sum(nvector);
Phi     =zeros(M,N);
temp    =0;

for i=1:length(nvector);
    ball                            =normc(randn(M,nvector(i)));
    d                               =diag(rand(nvector(i),1).^(1/M));
    ball                            =ball*d + c(:,i) * ones(1,length(d)) ;
    Phi(:,temp+1:temp+nvector(i))   =ball;
    temp                            =temp+nvector(i);
end

% create_dist_matrix--------------
D=zeros(N);

for i=1:N

    for j=1:N
        D(i,j)=norm(Phi(:,i)-Phi(:,j))^2;
    end
end

% step1--------------
nu=zeros(N,1);

for i=1:N
    nu(i)=norm(Phi(:,i))^2;
end

% step 4--------------
mu{1,k}=[];

for a=1:k
    na      =nvector(a);
    nua     =nu(n_sum(a)+1:n_sum(a+1));
    Phia    =Phi(:,n_sum(a)+1:n_sum(a+1));
    
    temp1   =step2_old(na,na,ones(na,1),Phia,Phia,nua,nua);
    temp2   =2/na*temp1;
    temp3   =ones(1,na)*temp1;
    temp3   =ones(na,1)*temp3/na^2;
    mua     =(temp3-temp2)/2;
    mu(a)   ={mua};
end


%multiple steps--------------
Mab1s{k,k}  =[];
z           =realmax;

for a=1:k

    for b=1:k
    
        if a~=b
        
            %step 5--------------
            mua         =mu{a};
            mub         =mu{b};
            na          =nvector(a);
            nb          =nvector(b);
            nua         =nu(n_sum(a)+1:n_sum(a+1));
            nub         =nu(n_sum(b)+1:n_sum(b+1));
            Phia        =Phi(:,n_sum(a)+1:n_sum(a+1));
            Phib        =Phi(:,n_sum(b)+1:n_sum(b+1));
            temp1       =mub'*ones(nb,1);
            temp2       =ones(1,nb)*ones(nb,1);
            temp3       =step2_old(na,nb,ones(nb,1),Phia,Phib,nua,nub);
            temp1       =ones(na,1)*temp1;
            temp2       =mua*temp2;
            Mab1        =temp3+temp2+temp1;
            Mab1s(a,b)  ={Mab1};
            
            %step 6--------------
            candidate   =2*nvector(a)*min(Mab1)/(nvector(a)+nvector(b));
            test        =(candidate<=z);
            z           =z*(1-test)+candidate*test;
        end
    end
end
clear mu mua mub Phia Phib candidate test;

%multiple steps--------------
uabs{k,k}=[];
rhoabs=zeros(N,1);

for a=1:k

    for b=1:k
    
        if a~=b
        
            %step 7--------------
            na          =nvector(a);
            nb          =nvector(b);
            uabs(a,b)   ={Mab1s{a,b}-z*(na+nb)/(2*na)*ones(na,1)};
            
            %step 8--------------
            temp        =uabs{a,b}'*ones(na,1);
            rhoabs(a,b) =temp;
            rhoabs(b,a) =temp;
        end
    end
end
clear Mab1s;

%step 9 (make all of B)
B=zeros(N);

for a=1:k

    for b=1:k
    
        if a~=b
            uab         =uabs{a,b};
            uba         =uabs{b,a};
            temp1       =uab*uba';
            rhoba       =rhoabs(b,a);
            temp_array  =temp1/rhoba;
            B(n_sum(a)+1:n_sum(a+1),n_sum(b)+1:n_sum(b+1))=temp_array;
            clear temp_array;
        end
    end
end
clear uabs rhoabs 

%step 10
Plambdaperp=eye(N);

for i=1:k
    temp                                                =zeros(N);
    temp(n_sum(i)+1:n_sum(i+1),n_sum(i)+1:n_sum(i+1))   =1/nvector(i);
    Plambdaperp                                         =Plambdaperp-temp;
end

eval            =max(real(eig(Plambdaperp*(B-D)*Plambdaperp)));
Theorem7check   =(eval<=z);

end

%--------------------------------------------------------------------------
%-------------------------------sub functions------------------------------
%--------------------------------------------------------------------------

%---------------

function output=step2_old(na,nb,x,Phia,Phib,nua,nub)

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

%---------------