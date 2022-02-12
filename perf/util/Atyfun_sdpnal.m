%%***********************************************************************
%% Atyfun: compute sum_{k=1}^m yk*Ak. 
%%
%%  Q = Atyfun(blk,At,y);
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

  function Q = Atyfun_sdpnal(blk,At,y)

  Q = cell(size(blk,1),1);
%%
  if iscell(At)
     for p = 1:size(blk,1)
        pblk = blk(p,:);
        n = pblk{1,2};
        if strcmp(pblk{1},'s')
           if (isempty(At{p}))
               Q{p} = sparse(n,n);
           else
               %Q{p} = mexsmat_sdpnal(pblk,At{p,1}*y);
               Q{p} = mexsmat(pblk,At{p,1}*y);
           end
        elseif strcmp(pblk{1},'l') || strcmp(pblk{1},'u')
           if (isempty(At{p}))
               Q{p} = zeros(n,1);
           else
               Q{p} = At{p,1}*y;
           end
        end 
     end
  else
     n = blk{1,2};
     if strcmp(blk{1,1},'s')
         if (isempty(At))
             Q = sparse(n,n);
         else
             %Q = mexsmat_sdpnal(blk,At*y);
             Q = mexsmat(blk,At*y);
         end
     else %%if strcmp(blk{1,1},'l') | strcmp(blk{1,1},'u')
         if (isempty(At))
             Q = sparse(n,n);
         else
             Q = At*y;
         end        
     end
  end
%%********************************************************* 

