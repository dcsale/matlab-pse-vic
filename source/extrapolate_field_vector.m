function uf = extrapolate_field_vector(SIM, MESH, uf)
% This routine extrapolates field values into the ghost layer

%% for now, just use nearest neighbor extrapolation.
% west (-x) and east (+x)
for i = 1:SIM.mbc
    uf{1}(i,:,:) = uf{1}(SIM.mbc+1,:,:);
    uf{2}(i,:,:) = uf{2}(SIM.mbc+1,:,:);
    uf{3}(i,:,:) = uf{3}(SIM.mbc+1,:,:);
end
for i = SIM.mbc+MESH.NX(1):MESH.NX(1)+2*SIM.mbc
    uf{1}(i,:,:) = uf{1}(SIM.mbc+MESH.NX(1),:,:);
    uf{2}(i,:,:) = uf{2}(SIM.mbc+MESH.NX(1),:,:);
    uf{3}(i,:,:) = uf{3}(SIM.mbc+MESH.NX(1),:,:);
end

% south (-y) and north (+y) 
for j = 1:SIM.mbc
    uf{1}(:,j,:) = uf{1}(:,SIM.mbc+1,:);
    uf{2}(:,j,:) = uf{2}(:,SIM.mbc+1,:);
    uf{3}(:,j,:) = uf{3}(:,SIM.mbc+1,:);
end
for j = SIM.mbc+MESH.NX(2):MESH.NX(2)+2*SIM.mbc
    uf{1}(:,j,:) = uf{1}(:,SIM.mbc+MESH.NX(2),:);
    uf{2}(:,j,:) = uf{2}(:,SIM.mbc+MESH.NX(2),:);
    uf{3}(:,j,:) = uf{3}(:,SIM.mbc+MESH.NX(2),:);
end

% bottom (-z) and top (+z)
for k = 1:SIM.mbc
    uf{1}(:,:,k) = uf{1}(:,:,SIM.mbc+1);
    uf{2}(:,:,k) = uf{2}(:,:,SIM.mbc+1);
    uf{3}(:,:,k) = uf{3}(:,:,SIM.mbc+1);
end
for k = SIM.mbc+MESH.NX(3):MESH.NX(3)+2*SIM.mbc
    uf{1}(:,:,k) = uf{1}(:,:,SIM.mbc+MESH.NX(3));
    uf{2}(:,:,k) = uf{2}(:,:,SIM.mbc+MESH.NX(3));
    uf{3}(:,:,k) = uf{3}(:,:,SIM.mbc+MESH.NX(3));
end

end % function


