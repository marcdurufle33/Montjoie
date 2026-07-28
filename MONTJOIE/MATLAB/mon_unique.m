 function [a,b]=mon_unique(sig,seuil)
 if(nargin==1)
     seuil=1e-6;
 end
%  
%  
% [a1,b1]=unique(sig,'last');
% c1=calcule_vitesse(a1',1);
% c1=c1-c1(1);
% 
% [m,ind]=find(abs(c1)<seuil);
% 

%% methode 1 
Delta = sig(2:end)-sig(1:end-1);
[index] = find(Delta<0);
if(isempty(index))
    a=sig;
    b=[1:length(sig)];
    return;
end
a=sig(1:index(1));
b=[1:index(1)];

for n=1:length(index)
   [iend] =  find(sig(index(n)+1:end)>a(end));
   fin = length(sig);
   try index(n+1);
       fin = index(n+1);
   catch 
       % rien
   end
   a=[a ;sig(index(n)+iend(1):fin)];
   b=[b index(n)+iend(1):fin];
end
%% methode 2
% compt=1;
% b=zeros(length(sig),1);
% a(1)=sig(1);
% for i=1:length(sig)-1
%     if(sig(i+1)>a(compt))
%         
%         compt=compt+1;
%         b(compt)=i+1;
%         a(compt)=sig(i+1);
%     end
% end