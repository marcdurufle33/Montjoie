function ErreurRelative = erreurMediane(varargin)

ratio = 0.9;
nb_vec = nargin/2;
norme_diff = 0;
norme_abs = 0;
for i = 1:nb_vec
  V1 = reshape(varargin{2*i-1}, numel(varargin{2*i-1}), 1);
  W1 = reshape(varargin{2*i}, numel(varargin{2*i}), 1);
  threshold = 1e-7*max(abs(W1));
  V1 = V1.*(abs(V1)>threshold).*(abs(W1)>threshold);
  W1 = W1.*(abs(V1)>threshold).*(abs(W1)>threshold);
  diff = sort(abs(V1-W1));
  val = sort(abs(V1));  
  offset = min(find(val > 0));
  N = numel(val) - offset;
  norme_diff = norme_diff + diff(floor(offset + ratio*N));
  norme_abs = norme_abs + val(floor(offset + ratio*N));
end

ErreurRelative = norme_diff/norme_abs;
