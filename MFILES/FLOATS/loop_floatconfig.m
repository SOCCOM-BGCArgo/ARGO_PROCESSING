
ff = [20804
20328
20532
20135
20592
20084
20912
20060];

for i = 1:length(ff)
    mbari_fn = ['ua',num2str(ff(i))];
    tf = build_apex_config(mbari_fn)
end
