% Author:       Takayuki Iguchi
% Filename:     figure3_right.m
% Last edited:  14 May 2016
% Description:  This function implements the procedure outlined in Figure
%               3(right) in [1]. Below is the description as listed in the
%               paper: 
% 
%               "For each Delta = 2 : 0.1 : 3 and k= 2 : 2 : 20, consider 
%               the unit balls in R^20 centered at 
%               \{ Delta * e_i / sqrt(2) \}_{i=1}^k where e_i denotes the 
%               ith identity basis element. Draw 100 points uniformly from 
%               each ball, and test if the resulting data points satisfy 
%               (12). After performing 10 trials of this experiment for 
%               each (Delta, k), we shade the corresponding pixel according
%               to the proportion of successful trials (white means every 
%               trial satisfied (12)). This plot indicates that our 
%               certificate (12) is to blame for Theorem 9’s dependence on 
%               k."
%               
%               See [1] for more details. 
%Inputs:
%               -d0: 
% 
%               A parameter used to determine the minimum value of Delta 
%               with  which we generate the points. For the paper, this 
%               value is 1.9 (which amounts to having a minimum value of 
%               Delta being 2 using an increment of 0.1)
% 
%               -m: 
% 
%               The dimension of the data to be generated. In the paper
%               this number was 20.
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
%               The number of trials for each combination of Delta and m
%               tested.
% 
%               -n: 
% 
%               Number of points to be drawn from each ball. 
%Outputs:
%   a phase transition described in the Program Description with the 
%   horizontal axis representing the minimum distance between clusters and 
%   the vertical axis representing the dimension m.

function figure3_right(d0,m,dcount,dsteps,trials,n)

data=zeros(m-1,dcount,trials);

for d=1:dcount,
    
    for k=2:m
    
        if mod(k,2)==0
            count=k-1;
        
            parfor ts=1:trials,
                data(count,d,ts)=kmeans_certificate_balls(d0+d*dsteps,n*ones(1,k),m)/trials;
            end
            disp([d k])
        end
    end
end
data    =data(1:2:m-1,:,:);
c0      =0;
ccount  =size(data,1);
cstep   =2;
data2   =zeros(ccount,dcount);

for d=1:dcount
    
    for increment_c=1:ccount
        data2(increment_c,d)=sum(data(increment_c,d,:));
    end
end
title_name=['Phase transition for fast dual cert in Relaxion, n',num2str(n),',',num2str(trials),' trials, dim=',num2str(m)];
make_phase_transition(d0,c0,m,dcount,ccount,cstep,dsteps,trials,data2,n,title_name)
end

%--------------------------------------------------------------------------
%----------------------------sub functions---------------------------------
%--------------------------------------------------------------------------

function make_phase_transition(d0,c0,m,dcount,ccount,cstep,dsteps,trials,data2,n,title_name)


success_prob=flipud(data2);
disp(success_prob);
clims=[0 1];
imagesc(success_prob,clims)
colormap(gray)


xticklabels = d0+ dsteps*[1:dcount];
xticks = linspace(1, size(data2, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)


yticklabels = (c0+ cstep*[1:ccount]);
yticks = linspace(1, size(data2, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))


title(title_name)

xlabel('\Delta')
ylabel('k')

end

