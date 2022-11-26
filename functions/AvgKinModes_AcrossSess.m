function [cdlateDel, cdMoveInit, cdResp] = AvgKinModes_AcrossSess(latents,taxis,meta)
tempCD.lateDel.R = []; tempCD.MoveInit.R = []; tempCD.Resp.R = [];
tempCD.lateDel.L = []; tempCD.MoveInit.L = []; tempCD.Resp.L = [];
for oo = 1:length(meta)
    tempCD.lateDel.R = [tempCD.lateDel.R, latents(oo).lateDel(:,1)];  tempCD.lateDel.L = [tempCD.lateDel.L, latents(oo).lateDel(:,2)];
    tempCD.MoveInit.R = [tempCD.MoveInit.R, latents(oo).MoveInit(:,1)];  tempCD.MoveInit.L = [tempCD.MoveInit.L, latents(oo).MoveInit(:,2)];
    tempCD.Resp.R = [tempCD.Resp.R, latents(oo).Resp(:,1)];  tempCD.Resp.L = [tempCD.Resp.L, latents(oo).Resp(:,2)];
end

cdlateDel.avg = NaN(length(taxis),2);  cdlateDel.std = NaN(length(taxis),2);
cdlateDel.avg(:,1) = mean(tempCD.lateDel.R,2,'omitnan'); cdlateDel.avg(:,2) = mean(tempCD.lateDel.L,2,'omitnan');
cdlateDel.std(:,1) = std(tempCD.lateDel.R,0,2,'omitnan'); cdlateDel.std(:,2) = std(tempCD.lateDel.L,0,2,'omitnan');

cdMoveInit.avg = NaN(length(taxis),2);  cdMoveInit.std = NaN(length(taxis),2);
cdMoveInit.avg(:,1) = mean(tempCD.MoveInit.R,2,'omitnan'); cdMoveInit.avg(:,2) = mean(tempCD.MoveInit.L,2,'omitnan');
cdMoveInit.std(:,1) = std(tempCD.MoveInit.R,0,2,'omitnan'); cdMoveInit.std(:,2) = std(tempCD.MoveInit.L,0,2,'omitnan');

cdResp.avg = NaN(length(taxis),2); cdResp.std = NaN(length(taxis),2);
cdResp.avg(:,1) = mean(tempCD.Resp.R,2,'omitnan'); cdResp.avg(:,2) = mean(tempCD.Resp.L,2,'omitnan');
cdResp.std(:,1) = std(tempCD.Resp.R,0,2,'omitnan'); cdResp.std(:,2) = std(tempCD.Resp.L,0,2,'omitnan');