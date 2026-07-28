function [x,y,line] = load1D_bb(fic_root)
  
fic = [fic_root '.bb'];
fid = fopen(fic,'r');
  n = fgetl(fid);
  line = fscanf(fid,'%f');
  fclose(fid);

  fic=[fic_root '.mesh'];
fid = fopen(fic,'r');
   n = fgetl(fid);
   n = fgetl(fid);
   n = fgetl(fid);
   n = fgetl(fid);
   n = str2num(fgetl(fid));
 [Mat] = fscanf(fid,'%f %f %d \n',[3,n]);
 
 fclose(fid);
x=Mat(1,:);
y=Mat(2,:);