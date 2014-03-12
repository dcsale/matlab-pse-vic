function Y = ReadParallelMatrix(X,file)
Y = X;   % Create output matrix.
myJ = global_ind(Y,2);  % Get local indices.
Yloc = local(Y);  % Create a local matrix.
for j=1:length(myJ)   % Loop over local indices.
    filej = [file '.' num2str(myJ(j)) '.mat'];  % Create filename.
    temp = load(filej,'Xj');  % Read data.
    Yloc(:,j) = temp.Xj;
end
Y = put_local(Y,Yloc);   % Copy back to distributed matrix.
