


% 1.首先将数据结果保存在alm_result/mat 和alm_result/res（梯度结果保存在alm_result/mat_grad 和 alm_result/res_grad）； mat太大，所以我清空了，你往里面放就行。

% 2 运行demo_scale_res.m得到scale后的结果，并保存在alm_result/re-scale_res (梯度结果保存在alm_result/re-scale_res_grad), 这一步很耗时

%3 运行print_table, 得到统计表和详细的表 （目前是保存了eta_K* 没有保存eta_d, 但是因为计算准则的时候，我投影了S，所以eta_K* = 0）， 有几个参数可以适当调节： app = 1.2 以及 判断成功的标准1e-5

%4. 运行plot_perf 得到profile图  nd = 3 全部画， nd =2 只画切换版本

