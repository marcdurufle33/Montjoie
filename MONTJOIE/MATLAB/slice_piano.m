function hout = slice_piano(xin,yin,zin,v,sx,sy,sz,xp,yp,zp,vp,zplaque)

% [cax,args,nargs] = axescheck(varargin{:});
% error(nargchk(4,8,nargs,'struct'));
% 
% cax = newplot(cax);
% hold_state = ishold(cax);
axis_slice = 1; % Slice only along axes?
method = 'linear';
cax=[];
cax = newplot(cax);

[ny,nx,nz] = size(v);
zp=zplaque*ones(size(xp));
if min(size(xin))==1 || min(size(yin))==1 ||  min(size(zin))==1,
    [x,y,z] = meshgrid(xin,yin,zin);
else
    [x,y,z] = deal(xin,yin,zin);
end
[ny,nx,nz] = size(v);
method = 'nearest';
hold_state = 0;



if axis_slice,
    h = [];
    [xi,yi,zi] = meshgrid(sx,y(:,1,1),z(1,1,:));
    
    vi = interp3(x,y,z,v,xi,yi,zi,method);
    for i = 1:length(sx)
          
        h = [h; surface( ...
            reshape(xi(:,i,:),[ny nz]),reshape(yi(:,i,:),[ny nz]), ...
            reshape(zi(:,i,:),[ny nz]),reshape(vi(:,i,:),[ny nz]), ...
            'parent',cax)];
    end

    [xi,yi,zi] = meshgrid(x(1,:,1),sy,z(1,1,:));
    vi = interp3(x,y,z,v,xi,yi,zi,method);
    for i = 1:length(sy)
        h = [h; surface( ...
            reshape(xi(i,:,:),[nx nz]),reshape(yi(i,:,:),[nx nz]),...
            reshape(zi(i,:,:),[nx nz]),reshape(vi(i,:,:),[nx nz]), ...
            'parent',cax)];
    end

%     [xi,yi,zi] = meshgrid(x(1,:,1),y(:,1,1),sz);
%     vi = interp3(x,y,z,v,xi,yi,zi,method);
%     for i = 1:length(sz)
%         h = [h; surface(xi(:,:,i),yi(:,:,i),zi(:,:,i),vi(:,:,i), ...
%             'parent',cax)];
%     end
    h=[h; surface(xp,yp,zp,vp,'parent',cax)];
else
    vi = interp3(x,y,z,v,sx,sy,sz,method);
    h = surf(sx,sy,sz,vi,'parent',cax);
end

if nargout > 0
    hout = h;
end
if ~hold_state
    view(cax,3), grid(cax,'on')
end

% Use ISFINITE to make sure no NaNs or Infs get passed to CAXIS
u=v(isfinite(v)); u = u(:);
caxis(cax,[min(u) max(u)])
% Signal to the world that we have created a new plot:
plotdoneevent(cax,h);

function str=id(str)
str = ['MATLAB:slice:' str];
