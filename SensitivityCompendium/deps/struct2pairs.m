function C=struct2pairs(S)
%Source: https://www.mathworks.com/matlabcentral/answers/469700-how-to-pass-a-struct-as-name-value-pairs-to-a-function#answer_381511
%Turns a scalar struct S into a cell of string-value pairs C
%
%  C=struct2pairs(S)
%
%If S is a cell already, it will be returned unchanged.
if iscell(S)
 C=S; return
elseif length(S)>1
    error 'Input must be a scalar struct or cell';
end
C=[fieldnames(S).'; struct2cell(S).'];
C=C(:).';