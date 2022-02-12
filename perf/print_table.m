close all;
clear all;

% compare scale result(alm_result/re-scale_res) with sdpnal result (sdpnal_result/res)

%res = 'res'; % non-scale-res
res = '../../scca_result/random/res2/'; %scale res
%res  = 'res-scale_res_1' % scale  res with normS





alm_res_name = dir(res); 

nid = length(alm_res_name)-2;



addpath 'util'



%data_path = 'sdp_data/ncm/';




table_str = '';table_str_s = '';
for i = 1:nid
    basename1 = alm_res_name(i+2).name;
    out = load([res,basename1]);
    
    
    [base,type] = file_split(basename1);
   
   % base = strrep(basename1,',mat','');
    r = out.r; p = out.p; q = out.q; lambda = out.lambda1;
    
    
    %trace = [out.lossu_trace,out.lossu_trace, mean(out.rho_trace)];
    Lasso = [out.lossu_lasso,out.lossu_lasso, mean(out.rho_lasso)];
    Init = [out.lossu_init,out.lossu_init, mean(out.rho_init)];
   

    
   

    if out.lossu_trace < min([Lasso(1),Init(1)])
        traceU = sprintf('\\textbf{%.3f}', out.lossu_trace);
    else
        traceU = sprintf('%.3f', out.lossu_trace);
    end
    
    if out.lossv_trace < min([Lasso(2),Init(2)])
        traceV = sprintf('\\textbf{%.3f}', out.lossv_trace);
    else
        traceV = sprintf('%.3f', out.lossv_trace);
    end
    
    if mean(out.rho_trace) > min([Lasso(3),Init(3)])
        trace_rho = sprintf('\\textbf{%.2f}', mean(out.rho_trace));
    else
        trace_rho = sprintf('%.2f', mean(out.rho_trace));
    end
    

   
    
    table_str = [table_str, sprintf('%d/%d/%d/%d/%.3f   & %s & %s & %s   & %.3f & %.3f & %.2f & %.3f & %.3f & %.2f ', ...
        type,p,q,r,lambda,traceU,traceV,trace_rho,Lasso,Init)];
    table_str = [table_str ' \\ \hline' newline];
    

   
    
    
end


% print detail
newtable  = table_str;

newtable = strrep(newtable,'e-0','e-');
newtable = strrep(newtable,'e+0','e+');
newtable = strrep(newtable,'e-1','e-1');
newtable = strrep(newtable,'e+1','e+1');

save_path = strcat('result/','scca_random','.txt');
fid = fopen(save_path,'w+');
fprintf(fid,'%s',newtable);







