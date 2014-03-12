function [] = save_dmat(d, outputFile)
%SAVE_DMAT Save a dmat to disk as a single file
%   SAVE_DMAT(D, OUTPUTFILE) Saves the dmat D to disk in the file OUTPUTFILE.
%       The elements in dmat D are saved in column-major order.
%
%   NOTE: ONLY SUPPORTS 2D DMATS WITH BLOCK DISTRIBUTION

DEBUG = 0;

% CREATE OUTPUT FILE

% Get number of processors and current rank
global pMATLAB;
my_rank = pMATLAB.my_rank;
Ncpus = pMATLAB.comm_size;

% Get datatype of data stored in d
datatype = class(d.local);

% Get size of datatype in bytes
switch (datatype)
   case {'uint8', 'int8'}
      precision = 1;

   case {'uint16', 'int16'}
      precision = 2;

   case {'uint32', 'int32', 'single'}
      precision = 4;

   case {'uint64', 'int64', 'double'}
      precision = 8;

   otherwise
     error([datatype ' is not supported.']);
end

% Allocate disk space for dmat
% Should be done by only leader processor
if (my_rank == 0)
   % Get index ranges for all processors
   [inds1 inds2] = global_ranges(d);

   create_output(outputFile, {inds1 inds2}, datatype);
end

% Last rank must unlock on the output file to ensure that rank 0
% will be able to start writing to output file
if (my_rank == (Ncpus-1))
   unlock_output(outputFile);
end


% WRITE TO OUTPUT FILE

% Wait until output file is unlocked before writing to it
synch_output(outputFile);

% Open outputFile for writing
fid = fopen(outputFile, 'r+');
if (fid == -1) error(['Error opening ' outputFile]); end

% Iterate over columns
d_local = local(d);
[ind1 ind2] = global_ind(d);
[dim1 dim2] = size(d);
for ii = 1:length(ind2)
   % Compute the single index for first index in col ii
   ind = sub2ind([dim1 dim2], ind1(1), ind2(ii));

   % Reposition file pointer to first element in row
   fileLoc = (ind - 1) * precision;
   s = fseek(fid, fileLoc, 'bof');
   if (s == -1) error(['Error accessing ' outputFile]); end

   if (DEBUG)
      disp(['(' num2str(ind1(1)) ', ' num2str(ind2(ii)) ') = ' ...
                num2str(ind) ' = ' ...
                num2str(fileLoc)] );
   end

   % Write data to output file
   count = fwrite(fid, d_local(:, ii), datatype);
end

s = fclose(fid);
if (s == -1) error(['Error closing ' outputFile]); end

% UNLOCK OUTPUT FILE
unlock_output(outputFile);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE_OUTPUT
function create_output(outputFile, inds, datatype)

   DEBUG = 0;

   % Get indices for each dimension
   inds1 = inds{1};
   inds2 = inds{2};

   fid = fopen(outputFile, 'w');
   if (fid == -1) error(['Error opening ' outputFile]); end

   % Iterate through each rank
   for ii = 1:size(inds1, 1)
      % Compute lengths of each dimension for current rank
      length1 = inds1(ii, 3) - inds1(ii, 2) + 1;
      length2 = inds2(ii, 3) - inds2(ii, 2) + 1;

      % Allocate an array with the same size as the current rank's local data
      temp = zeros(length1, length2);

      if (DEBUG)
         disp(['Size of d on rank ' num2str(ii) ' is (' num2str([length1 length2]) ')']);
      end

      % Write array to output file
      count = fwrite(fid, temp, datatype);
   end

   s = fclose(fid);
   if (s == -1) error(['Error closing ' outputFile]); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNCH_OUTPUT
function synch_output(outputFile)

   DEBUG = 0;

   % Get number of processors and current rank
   global pMATLAB;
   my_rank = pMATLAB.my_rank;

   % Create lock file name
   lock_file = fullfile('MatMPI',[outputFile '.lock.' num2str(my_rank)]);

   % Loop until lock file appears
   lock_fid = fopen(lock_file,'r');
   while (lock_fid == -1)
      lock_fid = fopen(lock_file,'r');

      pause(0.1);
   end
   fclose(lock_fid);

   if (DEBUG)
      disp(['Obtained lock on ' outputFile]);
   end

   % Delete lock file
   delete(lock_file);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNLOCK_OUTPUT
function unlock_output(outputFile)

   DEBUG = 0;

   % Get number of processors and current rank
   global pMATLAB;
   my_rank = pMATLAB.my_rank;
   Ncpus = pMATLAB.comm_size;

   % Compute next processor's rank
   next_rank = my_rank + 1;

   % If current processor has last rank, set next rank to 0
   if (my_rank == Ncpus-1)
      next_rank = 0;
   end

   % Create lock file name
   lock_file = fullfile('MatMPI',[outputFile '.lock.' num2str(next_rank)]);

   % Create lock file.
   fclose(fopen(lock_file,'w'));

   if (DEBUG)
      disp(['Released lock on ' outputFile]);
   end

