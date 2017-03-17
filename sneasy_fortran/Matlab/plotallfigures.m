function hout=plotallfigures (h, filename)
%PLOTALLFIGURES
%        Plots all figures up to handle h on h.ps

% Klaus Keller 8/10/96 klkeller@phoenix.edu

figure(1);
cmd = strcat(['print -dpsc2 ' filename]);
eval(cmd)

cmd = strcat(['print -append -dpsc2 ' filename]);
for j = 2:h
   figure(j)
   eval(cmd)
end
hout = 0 ;
return



