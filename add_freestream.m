function uf = add_freestream(uf, ENV)

uf{1} = uf{1} + ENV.velFree(1);
uf{2} = uf{2} + ENV.velFree(2);
uf{3} = uf{3} + ENV.velFree(3);

end

