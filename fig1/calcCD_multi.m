function cd = calcCD_multi(objs,params,cond,epoch)

mu = cell(numel(objs),1);
sd = cell(numel(objs),1);
for sessix = 1:numel(objs)
    [mu{sessix},sd{sessix}] = calcMu_SD(objs{sessix},params,cond,epoch,params.alignEvent,sessix);
end  

mu = cell2mat(mu);
sd = cell2mat(sd);

cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)

end