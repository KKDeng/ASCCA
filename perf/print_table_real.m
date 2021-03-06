close all;
clear all;

% compare scale result(alm_result/re-scale_res) with sdpnal result (sdpnal_result/res)

%res = 'res'; % non-scale-res
res = '../../scca_result/real/res_s/'; %scale res
%res  = 'res-scale_res_1' % scale  res with normS





alm_res_name = dir(res); 

nid = length(alm_res_name)-2;



addpath 'util'



%data_path = 'sdp_data/ncm/';




table_str1 = '';table_str2 = '';table_str3 = '';table_str4 = '';

for k = 1:23
for i = 1:nid
    basename1 = alm_res_name(i+2).name;
    out = load([res,basename1]);

   r = out.r; p = out.p; q = out.q; 
   A = load('../data/geneData.mat','A');
   
       if p == size(A.A{k,1},1)
           type = k;
       else
           continue;
       end
   
    
    
    
    trace = out.Trace; Lasso = out.Lasso; Pena = out.Pena;

    rho_trace_t = out.rho_trace_t; 
    rho_lasso_t = out.rho_lasso_t;
    rho_pena_t = out.rho_pena_t;
    
   
    if mean(rho_trace_t(:)) > max([mean(rho_lasso_t(:)),mean(rho_pena_t(:))])
        if r == 1
            trace_rho_t = sprintf('\\textbf{%.2f}', rho_trace_t(:));
        elseif r == 2
            trace_rho_t = sprintf('\\textbf{(%.2f,%.2f)}', rho_trace_t(:));
        elseif r == 3
            trace_rho_t = sprintf('\\textbf{(%.2f,%.2f,%.2f)}', rho_trace_t(:));
        else
            trace_rho_t = sprintf('\\textbf{(%.2f,%.2f,%.2f,%.2f)}', rho_trace_t(:));
        end
    else
        if r == 1
            trace_rho_t = sprintf('%.2f', rho_trace_t(:));
        elseif r == 2
            trace_rho_t = sprintf('(%.2f,%.2f)', rho_trace_t(:));
        elseif r == 3
            trace_rho_t = sprintf('(%.2f,%.2f,%.2f)', rho_trace_t(:));
        else
            trace_rho_t = sprintf('(%.2f,%.2f,%.2f,%.2f)', rho_trace_t(:));
        end
    end 
    
    
    if mean(trace(:)) > max([mean(Lasso(:)),mean(Pena(:))])
        if r == 1
            trace_rho = sprintf('\\textbf{%.2f}', trace(:));
        elseif r == 2
            trace_rho = sprintf('\\textbf{(%.2f,%.2f)}', trace(:));
        elseif r == 3
            trace_rho = sprintf('\\textbf{(%.2f,%.2f,%.2f)}', trace(:));
        else
            trace_rho = sprintf('\\textbf{(%.2f,%.2f,%.2f,%.2f)}', trace(:));
        end
    else
        if r == 1
            trace_rho = sprintf('%.2f', trace(:));
        elseif r == 2
            trace_rho = sprintf('(%.2f,%.2f)', trace(:));
        elseif r == 3
            trace_rho = sprintf('(%.2f,%.2f,%.2f)', trace(:));
        else
            trace_rho = sprintf('(%.2f,%.2f,%.2f,%.2f)', trace(:));
        end
    end
    
    
    
    
    
    
    if r == 1
        table_str1 = [table_str1 'chome',mat2str(type),'& '];
        table_str1 = [table_str1, sprintf('(%d, %d)   & %s & %s  & %.2f & %.2f & %.2f & %.2f ', ...
        p,q,trace_rho_t, trace_rho,rho_lasso_t,Lasso,rho_pena_t,Pena)];
    table_str1 = [table_str1 ' \\ \hline' newline];
    elseif r == 2
        table_str2 = [table_str2 'chome',mat2str(type),'& '];
    table_str2 = [table_str2, sprintf('(%d, %d)   & %s & %s  & (%.2f,%.2f) & (%.2f,%.2f) & (%.2f,%.2f) & (%.2f,%.2f)', ...
        p,q,trace_rho_t,trace_rho,rho_lasso_t,Lasso,rho_pena_t,Pena)];
    table_str2 = [table_str2 ' \\ \hline' newline];
    elseif r == 3
        table_str3 = [table_str3 'chome',mat2str(type),'& '];
        table_str3 = [table_str3, sprintf('(%d, %d)   & %s  & %.2f/%.2f/%.2f & %.2f/%.2f/%.2f ', ...
        p,q,trace_rho,Lasso,Pena)];
    table_str3 = [table_str3 ' \\ \hline' newline];
    else
        table_str4 = [table_str4 'chome',mat2str(type),'& '];
         table_str4 = [table_str4, sprintf('(%d, %d)   & %s  & %.2f/%.2f/%.2f/%.2f & %.2f/%.2f/%.2f/%.2f ', ...
        p,q,trace_rho,Lasso,Pena)];
    table_str4 = [table_str4 ' \\ \hline' newline];
    end
    
    

   
    
    
end
end
table_str  = [table_str1,table_str2,table_str3,table_str4];

% print detail
newtable  = table_str;

newtable = strrep(newtable,'e-0','e-');
newtable = strrep(newtable,'e+0','e+');
newtable = strrep(newtable,'e-1','e-1');
newtable = strrep(newtable,'e+1','e+1');

save_path = strcat('result/','scca_real_s','.txt');
fid = fopen(save_path,'w+');
fprintf(fid,'%s',newtable);







