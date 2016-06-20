% Author:       Takayuki Iguchi
% Filename:     figure3_left.m
% Last edited:  8 May 2016
% Description:  This script creates Figure 3 (left) in [1]. Below is the
%               caption for Figure 3 (left) in the paper:
%
%               "Take two unit disks in R^2 with centers on the x-axis at
%               distance 2.08 apart. Let x0 denote the smallest possible
%               x-coordinate in the disk on the right. For each disk, draw
%               N/2 = 50000 points uniformly at random from the
%               perimeter. Given theta, cluster the points according to whether
%               the x-coordinate is smaller than x0 + theta. When d = 0, this
%               clustering gives the planted clusters, and the k-means
%               objective (divided by N) is 1. We plot this normalized
%               k-means objective for theta in [0, 0.2]. Since N is large, this
%               curve is very close to its expected shape, and we see that
%               there are clusters whose k-means value is smaller than that
%               of the planted clustering."
%               
%               The script has a commented section which gives a
%               visualization of the clustered points for each value of theta.
% 
% Key Parameters:
%               -N:
%
%               Number of points to be clustered.
%           
%               -epsilon: 
%               
%               A parameter so that the unit disks have distance 2+epsilon
%               between disk centers.
% 
% Outputs:
%               A figure plotting the k-means objective for varying
%               inter-cluster center distances.
%
% Documentation:
%
% [1] Iguchi, Mixon, Peterson, Villar. Probably Certifiably Correct k-means
%       Clustering
% -------------------------------------------------------------------------

%Initialize important parameters
N           =100000;%<--User can edit this value
epsilon     =0.08;%<--User can edit this value
n           =[N/2,N/2];
c           =[1,0;-1-epsilon,0];
theta_min   =0;
theta_max   =0.2;
theta_step  =0.001;
plot_data   =zeros(size(theta_min:theta_step:theta_max,2),1);

%generate data
p               =zeros(N,2);
angles          =2*pi*rand(N/2,1);
u               =[cos(angles)+1,sin(angles)];
p(1:N/2,:)      =u;
angles          =2*pi*rand(N/2,1);
u               =[cos(angles)-1-epsilon,sin(angles)];
p(1+N/2:N,:)    =u;

%initialize for loop
count=1;

for theta=theta_min:theta_step:theta_max
    
    %pull points from A
    A   =p(p(:,1)<theta,:);
    B   =p(p(:,1)>=theta,:);
    
%     %visualize points being clustered
%     figure(1)
%     scatter(A(:,1),A(:,2))
%     hold on
%     scatter(B(:,1),B(:,2))
%     hold off
    
    %calulate emperical centers
    mu_A        =sum(A,1)/numel(A(:,1));
    mu_B        =sum(B,1)/numel(B(:,1));
    mu_A_array  =zeros(numel(A),2);
    
    %calculate f(theta)
    sq_diff_A=0;
    
    for j=1:numel(A(:,1))
        sq_diff_A=norm(A(j,:)-mu_A)^2+sq_diff_A;
    end
    sq_diff_B=0;
    
    for j=1:numel(B(:,1))
        sq_diff_B=norm(B(j,:)-mu_B)^2+sq_diff_B;
    end
    plot_data(count)=(sq_diff_A+sq_diff_B)/N;
    count=count+1;
end
figure(2);
plot(theta_min:theta_step:theta_max,plot_data,theta_min:theta_step:theta_max,ones(size(theta_min:theta_step:theta_max)))
title(['k-means objective at epsilon=',num2str(epsilon)])


