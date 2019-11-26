% JOSH - 11/12/17
%
% Select bootstrap index sorted without repeats

function [index]=bootstrapindex(ndata,nperm)

% n=length(data);
% index=randi(n,n,1);
index = sort(randperm(ndata,nperm));
