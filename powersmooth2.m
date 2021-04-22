function vecs=powersmooth2(vec,order,weight)
% BM Friedrich, 19.12.2014

%This is a modified version, which is easier to understand. For the original version please check:
% https://www.mathworks.com/matlabcentral/fileexchange/48799-powersmooth?s_tid=srchtitle

% cost function to be minimized
% cost = @(vec,vecs) sum( (vec-vecs).^2 ) + ...
%                     weight*sum( diff(vecs,order).^2 );
% rewriting cost function as quadratic form: vecs.A.vecs+vecs.b+c->min

%An Gong, last modified 2021.04.21.

N=length(vec);
vecs=nan(N,1);
isRow=isrow(vec);
if isRow
  vec = vec';
end
badIndex=isnan(vec);
vec(badIndex)=[];
n = length(vec);
D = spdiags(ones(n,1),0, n, n)-spdiags(ones(n-1,1),-1, n, n); % D*vec is the first-order forward difference of vec plus an extral vec(1)
Ek = spdiags([zeros(order, 1); ones(n-order, 1)], 0, n, n); %remove the first n points
Dk = Ek*(D^order); % Dk*vec is exactly the k-th order forward difference of vec
A = speye(n)+weight*(Dk')* Dk; %this factor can be derived by taking the derivate of the cost function
v = A\vec;
vecs(~badIndex)=v;
if isRow
    vecs=vecs';
end
