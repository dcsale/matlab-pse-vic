function uf = extrapolate_field_scalar(SIM, MESH, uf)
% This routine extrapolates field values into the ghost layer

%% for now, just use nearest neighbor extrapolation.
% west (-x) and east (+x)
for i = 1:SIM.mbc
    uf(i,:,:) = uf(SIM.mbc+1,:,:);
end
for i = SIM.mbc+MESH.NX(1):MESH.NX(1)+2*SIM.mbc
    uf(i,:,:) = uf(SIM.mbc+MESH.NX(1),:,:);
end

% south (-y) and north (+y) 
for j = 1:SIM.mbc
    uf(:,j,:) = uf(:,SIM.mbc+1,:);
end
for j = SIM.mbc+MESH.NX(2):MESH.NX(2)+2*SIM.mbc
    uf(:,j,:) = uf(:,SIM.mbc+MESH.NX(2),:);
end

% bottom (-z) and top (+z)
for k = 1:SIM.mbc
    uf(:,:,k) = uf(:,:,SIM.mbc+1);
end
for k = SIM.mbc+MESH.NX(3):MESH.NX(3)+2*SIM.mbc
    uf(:,:,k) = uf(:,:,SIM.mbc+MESH.NX(3));
end

end % function