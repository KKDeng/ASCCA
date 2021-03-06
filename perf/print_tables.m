close all;
clear all;

% compare scale result(alm_result/re-scale_res) with sdpnal result (sdpnal_result/res)

%res = 'res'; % non-scale-res
res = '../../scca_result/random/res_s/'; %scale res
%res  = 'res-scale_res_1' % scale  res with normS

types = 4; 



alm_res_name = dir(res); 

nid = length(alm_res_name)-2;



addpath 'util'



%data_path = 'sdp_data/ncm/';




table_str1 = ''; table_str2 = '';table_str3 = ''; table_str4 = ''; 
for i = 1:nid
    basename1 = alm_res_name(i+2).name;
    out = load([res,basename1]);
    
    
    [base,type] = file_split(basename1);
    if type ~= types 
        continue; 
    end
   % base = strrep(basename1,',mat','');
    r = out.r; p = out.p; q = out.q; 
    
    
    trace = out.Trace; Lasso = out.Lasso; Pena = out.Pena;

    
   

    if trace(1) < min([Lasso(1),Pena(1)])
        traceU = sprintf('\\textbf{%.3f}', trace(1));
    else
        traceU = sprintf('%.3f', trace(1));
    end
    
    if trace(2) < min([Lasso(2),Pena(2)])
        traceV = sprintf('\\textbf{%.3f}', trace(2));
    else
        traceV = sprintf('%.3f', trace(2));
    end
    
    if mean(trace(3:end)) > max([mean(Lasso(3:end)),mean(Pena(3:end))])
        if r == 1
            trace_rho = sprintf('\\textbf{%.2f}', trace(3:end));
        elseif r == 2
            trace_rho = sprintf('\\textbf{(%.2f,%.2f)}', trace(3:end));
        elseif r == 3
            trace_rho = sprintf('\\textbf{(%.2f,%.2f,%.2f)}', trace(3:end));
        else
            trace_rho = sprintf('\\textbf{(%.2f,%.2f,%.2f,%.2f)}', trace(3:end));
        end
    else
        if r == 1
            trace_rho = sprintf('%.2f', trace(3:end));
        elseif r == 2
            trace_rho = sprintf('(%.2f,%.2f)', trace(3:end));
        elseif r == 3
            trace_rho = sprintf('(%.2f,%.2f,%.2f)', trace(3:end));
        else
            trace_rho = sprintf('(%.2f,%.2f,%.2f,%.2f)', trace(3:end));
        end
    end
    
    
   
    if r == 1
    table_str1 = [table_str1, sprintf('(%d,%d)   & %s & %s & %s   & %.3f & %.3f & %.2f & %.3f & %.3f & %.2f ', ...
        p,q,traceU,traceV,trace_rho,Lasso,Pena)];
    table_str1 = [table_str1 ' \\ \hline' newline];
    elseif r == 2
    table_str2 = [table_str2, sprintf('(%d,%d)   & %s & %s & %s   & %.3f & %.3f & (%.2f,%.2f) & %.3f & %.3f & (%.2f,%.2f) ', ...
        p,q,traceU,traceV,trace_rho,Lasso,Pena)];
    table_str2 = [table_str2 ' \\ \hline' newline];
    elseif r == 3
        table_str3 = [table_str3, sprintf('(%d,%d)   & %s & %s & %s   & %.3f & %.3f & (%.2f,%.2f,%.2f) & %.3f & %.3f & (%.2f,%.2f,%.2f) ', ...
        p,q,traceU,traceV,trace_rho,Lasso,Pena)];
    table_str3 = [table_str3 ' \\ \hline' newline];
    else
         table_str4 = [table_str4, sprintf('(%d,%d)   & %s & %s & %s   & %.3f & %.3f & (%.2f,%.2f,%.2f,%.2f) & %.3f & %.3f & (%.2f,%.2f,%.2f,%.2f) ', ...
        p,q,traceU,traceV,trace_rho,Lasso,Pena)];
    table_str4 = [table_str4 ' \\ \hline' newline];
    end
    
    
    
    

   
    
    
end

table_str = '';
table_str = [table_str,' \multicolumn{10}{|c|}{  $r = 1 $}','\\ \hline', newline];
table_str = [table_str,table_str1];
table_str = [table_str,'\multicolumn{10}{|c|}{  $r = 2 $}','\\ \hline', newline];
table_str = [table_str,table_str2];
table_str = [table_str,'\multicolumn{10}{|c|}{  $r = 3 $}','\\ \hline', newline];
table_str = [table_str,table_str3];
table_str = [table_str,'\multicolumn{10}{|c|}{  $r = 4 $}','\\ \hline', newline];
newtable = [table_str,table_str4];
% print detail
%newtable  = [table_str2,table_str3,table_str4];

newtable = strrep(newtable,'e-0','e-');
newtable = strrep(newtable,'e+0','e+');
newtable = strrep(newtable,'e-1','e-1');
newtable = strrep(newtable,'e+1','e+1');

save_path = strcat('result/','scca_random_s','_',mat2str(types),'.txt');
fid = fopen(save_path,'w+');
fprintf(fid,'%s',newtable);







