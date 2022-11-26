function [movement,cdlate,cdearly,cdgo] = AvgCDandMode_AcrossSess(rez,taxis,meta)
tempmove.R = [];
tempmove.L = [];
tempCD.late.R = []; tempCD.early.R = []; tempCD.go.R = [];
tempCD.late.L = []; tempCD.early.L = []; tempCD.go.L = [];
for oo = 1:length(meta)
    tempmove.R = [tempmove.R, rez(oo).moveselect(:,1)];  tempmove.L = [tempmove.L, rez(oo).moveselect(:,2)];
    tempCD.late.R = [tempCD.late.R, rez(oo).cdLate_latent(:,1)];  tempCD.late.L = [tempCD.late.L, rez(oo).cdLate_latent(:,2)];
    tempCD.early.R = [tempCD.early.R, rez(oo).cdEarly_latent(:,1)];  tempCD.early.L = [tempCD.early.L, rez(oo).cdEarly_latent(:,2)];
    tempCD.go.R = [tempCD.go.R, rez(oo).cdGo_latent(:,1)];  tempCD.go.L = [tempCD.go.L, rez(oo).cdGo_latent(:,2)];
end
movement.avg = NaN(length(taxis),2);  movement.std = NaN(length(taxis),2);
movement.avg(:,1) = mean(tempmove.R,2,'omitnan'); movement.avg(:,2) = mean(tempmove.L,2,'omitnan');
movement.std(:,1) = std(tempmove.R,0,2,'omitnan'); movement.std(:,2) = std(tempmove.L,0,2,'omitnan');

cdlate.avg = NaN(length(taxis),2);  cdlate.std = NaN(length(taxis),2);
cdlate.avg(:,1) = mean(tempCD.late.R,2,'omitnan'); cdlate.avg(:,2) = mean(tempCD.late.L,2,'omitnan');
cdlate.std(:,1) = std(tempCD.late.R,0,2,'omitnan'); cdlate.std(:,2) = std(tempCD.late.L,0,2,'omitnan');

cdearly.avg = NaN(length(taxis),2);  cdearly.std = NaN(length(taxis),2);
cdearly.avg(:,1) = mean(tempCD.early.R,2,'omitnan'); cdearly.avg(:,2) = mean(tempCD.early.L,2,'omitnan');
cdearly.std(:,1) = std(tempCD.early.R,0,2,'omitnan'); cdearly.std(:,2) = std(tempCD.early.L,0,2,'omitnan');

cdgo.avg = NaN(length(taxis),2); cdgo.std = NaN(length(taxis),2);
cdgo.avg(:,1) = mean(tempCD.go.R,2,'omitnan'); cdgo.avg(:,2) = mean(tempCD.go.L,2,'omitnan');
cdgo.std(:,1) = std(tempCD.go.R,0,2,'omitnan'); cdgo.std(:,2) = std(tempCD.go.L,0,2,'omitnan');