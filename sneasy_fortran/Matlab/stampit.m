function ht = stampit(h,shft,fn)
% PLTSTMP Adds a timestamp and additonal comment to a plot
%
% PLTSTMP adds a timestamp to the plot. If PLTSTMP is issued in a m-file,
% the full path and filename is also added in the footer of the plot.
%
% PLTSTMP(FID) adds the timestamp to figure FID instead of the current
% figure.
%
% The timestamp sometimes interferes with the x-axis label. To avoid
% this, PLTSTMP(FID,SHFT) moves child axes of FID vertically by SHFT
% (in normalized coordinates) before adding the timestamp.
%
% PLTSTMP(FID,SHFT,STR) adds STR to the footer instead of the m-file
% name. 
%
% H = PLTSTMP also returns the handle of the text string.
%
% Example:
% plot(pascal(5)), pltstmp

% aha, 3-jul-2004, fixed zoom problem and documentation
% aha, 15-oct-2003

if nargin<1, h=gcf; end

% shift original axes upward by shft
if nargin>1
  ax = get(gcf,'child');
  ax = findobj(gcf,'Type','Axes');
  for k = 1:length(ax)
    pp = get(ax(k),'Position');
    set(ax(k),'Position',pp+[0 shft 0 0])
  end
end

if nargin<3
  % no string given, retrieve m-file name
  [St,I]=dbstack;
  if length(St)>I
		fn = which(St(I+1,1).name);
    fn = [fn ' [' datestr(now) ']'];
  else
    fn = [' [' datestr(now) ']'];
  end
  in = 'non';
else
  % add timestamp to given string
  fn = [fn ' [' datestr(now) ']'];
  in = 'tex';
end

% delete a previous existing timestamp
set(0,'ShowHiddenHandles','On');
delete(findobj(h,'Tag','PLTSTMP-AXES'));
set(0,'ShowHiddenHandles','Off');

% create axes for text string and add string to the bottom of the
% figure
a = axes('parent',h,...
  'posi',[-0 -0 1 1],...
  'ydir','norm','xdir','norm',...
  'visible','off',...
  'Handlevisibility','off',...
  'Tag','PLTSTMP-AXES');
ht = text(.5,0,fn,'hor','center','vert','bot','fontsize',8,...
  'interpre',in,'parent',a);

% set(h,'currentaxes',oax)
if nargout==0, clear, end

