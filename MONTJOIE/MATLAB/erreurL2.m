function ErreurRelative = erreurL2(varargin)

epsil = 0
if (mod(nargin, 2) == 1)
  epsil = varargin(nargin)
end

nb_vec = nargin/2;
norme_diff = 0;
norme_abs = 0;
for i = 1:nb_vec
  V1 = reshape(varargin{2*i-1}, numel(varargin{2*i-1}), 1);
  W1 = reshape(varargin{2*i}, numel(varargin{2*i}), 1);
  threshold = epsil*max(abs(W1));
  V1 = V1.*(abs(V1)>threshold).*(abs(W1)>threshold);
  W1 = W1.*(abs(V1)>threshold).*(abs(W1)>threshold);
  norme_diff = norme_diff + sum((V1-W1).*conj(V1-W1));
  norme_abs = norme_abs + sum(W1.*conj(W1));
end

ErreurRelative = sqrt(norme_diff/norme_abs);
