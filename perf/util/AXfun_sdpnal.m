%%***********************************************************************
%% AXfun: compute AX(k) = <Ak,X>, k = 1:m
%%
%% AX = AXfun(blk,At,X);
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

  function AX = AXfun_sdpnal(blk,At,X)

  if iscell(At)
     m = size(At{1},2);  
     AX = zeros(m,1);    
     for p = 1:size(blk,1)
        pblk = blk(p,:);
        if strcmp(pblk{1},'s')
           if (isempty(At{p}))
              AX = [];
              AXtmp = [];
           elseif (length(pblk{2}) == 1)
              %AXtmp = (mexsvec_sdpnal(pblk,X{p})'*At{p,1})';
              AXtmp = (mexsvec(pblk,X{p})'*At{p,1})';
	         else
              %AXtmp = (mexsvec_sdpnal(pblk,sparse(X{p}))'*At{p,1})';
              AXtmp = (mexsvec(pblk,sparse(X{p}))'*At{p,1})';
           end
        elseif strcmp(pblk{1},'l') || strcmp(pblk{1},'u')
           if (isempty(At{p}))
              AX = [];
              AXtmp = [];
           else
              AXtmp = (X{p}'*At{p,1})';
           end
        end
	   AX = AX + AXtmp; 
     end
  else
     if strcmp(blk{1,1},'s')
        if (isempty(At))
           AX =[];
        elseif (length(blk{1,2})==1)
           %AX = (mexsvec_sdpnal(blk,X)'*At)';
           AX = (mexsvec(blk,X)'*At)';
	      else
           %AX = (mexsvec_sdpnal(blk,sparse(X))'*At)';
           AX = (mexsvec(blk,sparse(X))'*At)';
        end
     elseif strcmp(blk{1,1},'l') || strcmp(blk{1,1},'u')
        if (isempty(At))
           AX =[];
        else
           AX = (X'*At)';
        end
     end
  end
%%*********************************************************
