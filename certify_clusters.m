% Author:       Iguchi, Mixon, Peterson, Villar
% Filename:     certify_clusters.m
% Last edited:  14 May 2016
% Description:  This function provides a certificate of optimality for a 
%               given set of clusters using the procedure outlined
%               in Section 3.2 of [1] to test the sufficient condition for
%               k-means optimality in Theorem 7 in [1]. See [1] for more 
%               details. 
% Inputs:       
%               -Phi:
%               
%               A m x n array with the coordinates of n points in dimension
%               R^m. This array correspond to the points which clustering 
%               the algorithm certifies. They need to be ordered by
%               cluster.
% 
%               -nvector:
% 
%               A 1 x k array where the ith entry corresponds to the number
%               of points in the ith ball. The order must coincide with
%               the order of the points in Phi. The sum of the entries of
%               nvector must coincide with the number of columns of Phi.
% 
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

function [Theorem7check]=certify_clusters(Phi,nvector)

%initialize--------------
[m,n]   =size(Phi);
k       =max(size(nvector));
n_sum   =nvector(1)*ones(1,k);

for i=2:k
    n_sum(i)=nvector(i)+n_sum(i-1);
end
n_sum   =horzcat(0,n_sum);
N       =sum(nvector);
if n~=N
    disp('Error: the number of columns in Phi must coincide with the sum of the entries of nvector');
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
        end
    end
end
clear uabs rhoabs temp_array;

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