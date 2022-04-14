function ternlabel(atxt,btxt,ctxt);
% TERNLABEL - Ternary labels.
%
% ternlabel(atxt,btxt,ctxt)

w = 0.5;
h = 0.5/tan(pi/6);
d = 0.02;

text(0,h+d,atxt,'FontSize',12,'HorizontalAlignment','center','FontWeight','bold','VerticalAlignment','bottom');
text(-w-d,0,btxt,'FontSize',12,'HorizontalAlignment','right','FontWeight','bold');
text(w+d,0,ctxt,'FontSize',12,'HorizontalAlignment','left','FontWeight','bold');

return
