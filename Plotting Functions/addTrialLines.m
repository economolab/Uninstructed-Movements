function addTrialLines(col)
% Add lines at trial times
go = 0;
delstart = -0.9;
sampstart = delstart-1.3;
xline(go,'Color',col,'LineStyle','--')
xline(delstart,'Color',col,'LineStyle','--')
xline(sampstart,'Color',col,'LineStyle','--')
end
