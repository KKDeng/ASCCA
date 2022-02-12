function [name1,type] = file_split(basename)


if contains(basename,'rand')
    name1 = 'rand';
else
    name1 = 'real';
end
if  contains(basename,'rand_1_')
    type  =1;
elseif contains(basename,'rand_2_')
    type  =2;
elseif contains(basename,'rand_3_')
    type  =3;
elseif contains(basename,'rand_4_')
    type  =4;
end

% 
% if ~isempty(name1)
%     name1 = [name1 ' & '];
% end














