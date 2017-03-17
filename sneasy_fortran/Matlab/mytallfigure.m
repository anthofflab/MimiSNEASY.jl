function foo=mytallfigure(number,width)
% function myfigure(number)
% /home/klaus/wrk/matlab/mystuff
% Fri Jun 22 13:54:07 EDT 2001
% update: Mon Oct 29 13:31:59 EST 2001   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(number);
clf;
set(gcf,'PaperPosition',[1.75 1.5 5 8]); 
set(gcf,'DefaultLineLineWidth',width)
wysiwyg;
foo=number;



