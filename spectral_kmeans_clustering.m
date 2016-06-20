%Filename:      spectral_kmeans_clustering.m
%Author:        Takayuki Iguchi
%Last edit:     8 May 2016
%Description:   This program implements Algorithm 2 in [1].
%Inputs:
%               -Phi: 
%     
%               an m x N array of points to be clustered (in 2 clusters) 
%               where N denotes the number of points and m denotes the 
%               dimension of the space the points are in.
%               
%               -tol: 
% 
%               a positive number defining acceptable 2-norm error in 
%               finding the leading eigenvector of a matrix using the 
%               power method.
% 
%Outputs:
%               -IDX: 
% 
%               a N x 1 array of cluster labels for the N points in Phi.
% 
%Documentation:
% [1] Iguchi, Mixon, Peterson, Villar. Probably certifiably correct k-means
%       clustering
% -------------------------------------------------------------------------

function IDX=spectral_kmeans_clustering(Phi,tol)

[m,N]=size(Phi);

%subtract centroid from each column of Phi to produce Phi_0
centroid        =sum(Phi,2)/N;
centroid_array  =repmat(centroid,1,N);
Phi_0           =Phi-centroid_array;

%Compute the leading eigenvector y of Phi_0^T*Phi_0
x       =rand(m,1);
x       =x/norm(x,2);
cont    =true;
count   =1;

while cont
    xtemp   =Phi_0'*x;
    xnew    =Phi_0*xtemp;
    xnew    =xnew/norm(xnew);
    error   =norm(xnew-x);
    x       =xnew;
    
    if error<tol
        cont=false;
        
    elseif count>10000%some large number
        fprintf('too many interations\n')
        cont=false;
        
    else
        count=count+1;
        
    end
end
y=Phi_0'*x;

%Find theta that minimizes the kmeans objective of
%  ({i:y_i<theta},{i:y_i>=theta})

    %Sort the entries y_1<=...<=y_N in O(N log N) operations
    [~, idx_mapping]    =sort(y);
    sorted_Phi_0        =Phi_0(:,idx_mapping);

    %Iteratively compute s_1(i),s_1^c(i),s_2(i),s_2^c(i) for i=1:N-1 in
    %   O(mN) operations
    s_1         =zeros(m,N-1);
    s_1c        =zeros(m,N-1);
    s_2         =zeros(1,N-1);
    s_2c        =zeros(1,N-1);
    list_norms  =zeros(N,1);

    for i=1:N
        list_norms(i)=norm(sorted_Phi_0(:,i),2)^2;
    end
    s_1(:,1)    =sorted_Phi_0(:,1);
    s_2(1)      =list_norms(1);
    s_1c(:,N-1) =sorted_Phi_0(:,N);
    s_2c(N-1)   =list_norms(N);

    for i=2:N-1
        s_1(:,i)    =s_1(:,i-1)+sorted_Phi_0(:,i);
        s_2(i)      =s_2(i-1)+list_norms(i);
    end

    for i=N-2:-1:1
        s_1c(:,i)   =s_1c(:,i+1)+sorted_Phi_0(:,i);
        s_2c(i)     =s_2c(i+1)+list_norms(i);
    end

    %Compute v_1=0 and v_{i+1} for every i=1:N-2 in O(mN) operations
    v=zeros(N-1,1);

    for i=1:N-2
        v(i+1) = v(i)+2*s_2(i)-4*sorted_Phi_0(:,i+1)'*s_1(:,i)+2*i*norm(sorted_Phi_0(:,i+1),2)^2;
    end

    %Compute v_N^c=0 and v_^c{i-1} for every i=N:-1:2 in O(mN) operations
    vc=zeros(N-1,1);

    for i=N-1:-1:2
        vc(i-1) = vc(i)+2*s_2c(i)-4*sorted_Phi_0(:,i)'*s_1c(:,i)+2*(N-i)*norm(sorted_Phi_0(:,i),2)^2;
    end

    %Compute f(i) for every i=1:N-1 in O(N) operations
    f       =zeros(N-1,1);
    min_i   =1;
    min_f   =v(1)/1+vc(1)/(N-1);

    for i=2:N-1
        f(i)=v(i)/i+vc(i)/(N-i);

        if f(i)<=min_f
            min_i=i;
            min_f=f(i);
        end
    end

    %Find i that minimizes f(i) and output {i,...,i} and {i+1,...,N} in
    %   O(N) operations (see subsection above)

%Create IDX according to (A,B)<-({i:y_i<theta},{i:y_i>=theta})
IDX=(idx_mapping<=min_i)+1;
end
