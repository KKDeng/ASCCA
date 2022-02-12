%%plot_perf

close all;
clear all;
 name1 = 'random/';
 data_dir = ['../../',name1];
 alm_name = dir(data_dir);

%-------------------------------------------------------------------------

figpath =  'result/perf/';
name = strrep(name1,'/','');
res = '../../scca_result/random/res_s/'; 


%% 2018-07-version 

solver_range = {'Real','ASCCA','CoLaR','SCCA-PLS'};
colname_range = {'Real','ASCCA','CoLaR','SCCA-PLS'};
lines   = {'k-' 'b--' 'r-.' 'c:' 'g-*' 'm-' 'r-'};




re = {};       
    re{1} = load([res,'rand_1_300_300_1_0.8.mat']);
    re{2} = load([res,'rand_2_300_300_1_0.8.mat']);
    re{3} = load([res,'rand_3_300_300_1_0.8.mat']);
    re{4} = load([res,'rand_4_300_300_1_0.8.mat']);
    
    out = zeros(300,4);
for fidx = 1: 4

    
    
    fig = figure(fidx);
    
    out_u = [-re{fidx}.u_real,re{fidx}.u_trace,re{fidx}.u_lasso,re{fidx}.u_pena];
    
    out_v = [-re{fidx}.v_real,re{fidx}.v_trace,re{fidx}.v_lasso,re{fidx}.v_pena];
    
    

    x = [1:1:300];
    
    for s = 1:4
        
        plot(x,out_u(:,s),lines{s},'LineWidth',2,'MarkerSize',3);
        hold on
    end
    
   
    legend( colname_range, 'location','southeast');
    xlabel('not more than 2^x times worse than the best');
    ylabel('ratio of problems');
    title(sprintf('%s',perf_fig_name{fidx}));
    set(gca,'fontsize',15)
    print(fig , '-depsc2', strcat(figpath,name, '_perf_cmp_',perf_fig{fidx},'.eps'));
            
    
end

return
