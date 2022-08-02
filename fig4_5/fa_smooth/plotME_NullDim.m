function plotME_NullDim(obj,me,rez,sessix,nulldim,data_type)

%% plot motion energy and a potent space projection
potenttemp = rez(sessix).N_null;

metemp = me(sessix);


trix = 1:obj(sessix).bp.Ntrials;

medata = metemp.data(:,trix);

nulldata = potenttemp(:,trix,nulldim);

medata = mySmooth(medata(:),31);

nulldata = mySmooth(nulldata(:),31);

newtime = (1:numel(nulldata))./200;
f = figure; hold on;
f.Position = [355         328        1208         271];
% patchline(newtime+.120,medata,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
patchline(newtime,medata,'EdgeColor','k','EdgeAlpha',0.35,'LineWidth',2);
ix = medata > metemp.moveThresh;
z = medata;
z(~ix) = nan;
% plot(newtime+.120,z,'r','LineWidth',2)
plot(newtime,z,'r','LineWidth',2)


if strcmpi(data_type,'fa')
    plot(newtime,nulldata.*35 + 40,'b','LineWidth',1);
else
    plot(newtime,nulldata ./ 2 + 40 ,'b','LineWidth',1);
end




end