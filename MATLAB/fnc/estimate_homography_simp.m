function H = estimate_homography(X1, X2)
% function H = estimate_homography(X1, X2)
%
% Inputs:
% X1 = 2xN matrix of (x,y) locations in image 1
% X2 = 2xN matrix of (x,y) locations in image 2
%
% Outputs:
% H = 3x3 matrix of homography parameters
%
% Requires VLFeat library to already be in MATLAB path:
% run('/path/to/library/vlfeat-0.9.17/toolbox/vl_setup.m');

% --------------------------------------------------------------------
%                                         RANSAC with homography model
% --------------------------------------------------------------------

clear H score ok ;
numMatches = size(X1,2);
X1 = [X1; ones(1,size(X1,2))];
X2 = [X2; ones(1,size(X2,2))];
%nt = 100;
nt = 1; 
for t = 1:nt
  % estimate homography
  %subset = vl_colsubset(1:numMatches, 4) ;
  subset = 1:4;
  A = [] ;
  for i = subset
    A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
  end
  [U, S, V] = svd(A) ;
  H{t} = reshape(V(:,9),3,3) ;

  % score homography
  X2_ = H{t} * X1 ;
  du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
  dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
  ok{t} = (du.*du + dv.*dv) < 6*6 ;
  score(t) = sum(ok{t}) ;
end

[score, best] = max(score) ;
H = H{best} ;
ok = ok{best} ;

% --------------------------------------------------------------------
%                                                  Optional refinement
% --------------------------------------------------------------------
function err = residual(H)
 u = H(1) * X1(1,ok) + H(4) * X1(2,ok) + H(7) ;
 v = H(2) * X1(1,ok) + H(5) * X1(2,ok) + H(8) ;
 d = H(3) * X1(1,ok) + H(6) * X1(2,ok) + 1 ;
 du = X2(1,ok) - u ./ d ;
 dv = X2(2,ok) - v ./ d ;
 err = sum(du.*du + dv.*dv) ;
end

H = H / H(3,3) ;
opts = optimset('Display', 'none', 'TolFun', 1e-8, 'TolX', 1e-8) ;
H(1:8) = fminsearch(@residual, H(1:8)', opts) ;


end

