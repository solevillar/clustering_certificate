% Author:       Takayuki Iguchi
% Filename:     figure3_center.m
% Last edited:  14 May 2016
% Description:  This function implements the procedure outlined in Figure
%               3(center) in [1]. Below is the description as listed in the
%               paper: 
% 
%               "Take two intervals of width 2 in R, and let Delta denote  
%               the distance between the midpoints of these intervals. For 
%               each interval, draw 10 points at random from its endpoints,
%               and then run the k-means SDP. For each Delta = 2 : 0.1 : 5, 
%               after running 2000 trials of this experiment, we plot the 
%               proportion of trials for which the SDP relaxation was tight
%               (dashed line) and the proportion of trials for which the 
%               SDP recovered the planted clusters (solid line). In this 
%               case, cluster recovery appears to exhibit a phase 
%               transition at Delta = 4." 
%                 
%               See [1] for more details. 
% Inputs:       
%               -N: 
% 
%               The number of points to be clustered with N/2 points per
%               planted cluster. In the paper, this value was set to 20. 
% 
%               -d0: 
% 
%               A parameter used to determine the minimum value of Delta w/ 
%               which we generate the points. For the paper, this value is 
%               1.9 (which amounts to having a minimum value of Delta being 2 
%               using an increment of 0.1)
% 
%               -dcount: 
% 
%               The number of different values of Delta tested.
% 
%               -dsteps: 
% 
%               The increment used for varying Delta.
% 
%               -trials: 
% 
%               The number of trials for each value of Delta tested.
% Outputs:
%               -A figure comparing SDP integrality and SDP planted cluster
%               recovery for varying inter-cluster distances under the
%               distribution of points listed in the program description.
% 
%Documentation:
% 
% [1] Iguchi, Mixon, Peterson, Villar. Probably certifiably correct k-means
%       clustering
% -------------------------------------------------------------------------

function figure3_center(N,d0,dcount,dsteps,trials)

K                       =2;
test_vals_planted_array =zeros(dcount,trials);
test_vals_tight_array   =zeros(dcount,trials);

for d=1:dcount

    parfor ts=1:trials
    
        % generate data
        Delta   =d0+d*dsteps;
        p       =[2*binornd(1,0.5,N/2,1)-1;2*binornd(1,0.5,N/2,1)+1+(Delta-2)];

        % solve the SDP on the data using CVX
        sdp_clustering = kmeans_sdp_cvx(p,K);
        
        %check if the output is tight
        temp                            =validate_tight(sdp_clustering,10^-6);
        test_vals_tight_array(d,ts)     =temp/trials;
        temp2                           =validate_planted(sdp_clustering,p,10^-6);
        test_vals_planted_array(d,ts)   =temp2/trials;
        disp([d,ts])
    end
end
test_vals_planted   =sum(test_vals_planted_array,2);
test_vals_tight     =sum(test_vals_tight_array,2);

%plot the results
figure(11)
plot(d0+dsteps:dsteps:d0+dcount*dsteps,test_vals_tight,'--')
hold on
plot(d0+dsteps:dsteps:d0+dcount*dsteps,test_vals_planted,'black')
hold off
end

%--------------------------------------------------------------------------
%----------------------sub functions---------------------------------------
function testval=validate_planted(matrix,p,tol)

N               =size(p,1);
onea            =zeros(N,1);
onea(1:N/2)     =1;
oneb            =zeros(N,1);
oneb(N/2+1:N)   =1;

if norm(matrix*onea-onea)<tol && norm(matrix*oneb-oneb)<tol
    testval=1;
    
else
    testval=0;
end
end

function test_val=validate_tight(matrix,tol)

test_val    =0;
eigenvalues =eig(matrix);
eigenvalues =sort(eigenvalues);

if abs(eigenvalues(end-1)-1)<tol && abs(eigenvalues(end)-1)<tol

    if max(abs(eigenvalues(1:end-2)))<tol
        test_val=1;
    end
end
end