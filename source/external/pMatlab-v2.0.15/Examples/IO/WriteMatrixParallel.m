function WriteParallelMatrix(X,filebase)
myJ = global_ind(X,2);  % Get local indices.
Xloc = local(X);  % Get local matrix.
for j=1:length(myJ)  % Loop over local indices.
  Xj = Xloc(:,j);  % Get a vector.
  filej = [filebase '.' num2str(myJ(j)) '.mat'];  % Create filename.
  save(filej,'Xj');  % Save vector to a file.
end
