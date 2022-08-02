function plotME_PotentDim(obj,me,rez,sessix,potentdim,data_type)

%% plot motion energy and a potent space projection
potenttemp = rez(sessix).N_potent;

metemp = me(sessix);

trix = 1:obj(sessix).bp.Ntrials;

medata = metemp.data(:,trix);

potentdata = potenttemp(:,trix,potentdim);

medata = mySmooth(medata(:),31);

potentdata = mySmooth(potentdata(:),31);

newtime = (1:numel(potentdata))./200;
f = figure; hold on;
f.Position = [358         687        1198         294];
% patchline(newtime+.120,medata,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
patchline(newtime,medata,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
ix = medata > metemp.moveThresh;
z = medata;
z(~ix) = nan;
% plot(newtime+.120,z,'r','LineWidth',2)
plot(newtime,z,'r','LineWidth',2)


if strcmpi(data_type,'fa')
    plot(newtime,potentdata.*35 + 40,'b','LineWidth',1);
else
    plot(newtime,potentdata  + 40 ,'b','LineWidth',1);
end




end