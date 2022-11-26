function addTrialLines(col,met,obj)
% Add lines at trial times
go = 0;
trix = met.trialid{1}(1);
del = obj.bp.ev.goCue(trix)-obj.bp.ev.delay(trix);
delstart = 0-del;
sampstart = delstart-1.3;
xline(go,'Color',col,'LineStyle','--')
xline(delstart,'Color',col,'LineStyle','--')
xline(sampstart,'Color',col,'LineStyle','--')
end
