function [EToE,EToF]= tiConnect2D(EToV)
% function [EToE,EToF]= tiConnect2D(EToV)
% Purpose: quadrilateral face connect algorithm due to Toby Isaac
%
% Modification of routine so that it can be used for any type of elements
%
% by Allan P. Engsig-Karup, apek@imm.dtu.dk.
%
Nfaces = size(EToV,2);  % Number of faces of an element
K      = size(EToV,1);  % Number of elements
Nnodes = max(max(EToV));    % Number of vertices

% create list of all faces 1, then faces 2, ...
VNUM = zeros(Nfaces,2); % List of local face to local vertex connections
VNUM(:,1) = 1:Nfaces; VNUM(:,2) = [2:Nfaces 1];
fnodes = reshape(EToV(:,[VNUM(:,1) VNUM(:,2)]),Nfaces*K,2);
fnodes = sort(fnodes,2)-1;

% set up default element to element and Element to faces connectivity
EToE = (1:K)'*ones(1,Nfaces); EToF = ones(K,1)*(1:Nfaces);

% uniquely number each set of faces by their node numbers
id = fnodes(:,1)*Nnodes + fnodes(:,2)+1;
spNodeToNode=[id, (1:Nfaces*K)', EToE(:), EToF(:)];

% Now we sort by global face number.
sorted=sortrows(spNodeToNode,1);

% find matches in the sorted face list
[indices,dummy]=find( sorted(1:(end-1),1)==sorted(2:end,1) );

% make links reflexive
matchL = [sorted(indices,:)   ;sorted(indices+1,:)];
matchR = [sorted(indices+1,:) ;sorted(indices,:)];
% insert matches
EToE(matchL(:,2)) = matchR(:,3);
EToF(matchL(:,2)) = matchR(:,4);
return;
