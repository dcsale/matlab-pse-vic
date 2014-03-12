%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script implements a basic image convolution
% across multiple processors.
% To run, start Matlab and type:
%
%   eval( MPI_Run('blurimage',2,{}) );
%
% Or, to run a different machine type:
%
%   eval( MPI_Run('blurimage',2,{'machine1' 'machine2'}) );
%
% Output will be piped into to
%
%   MatMPI/blurimage.0.out
%   MatMPI/blurimage.1.out
%   ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatlabMPI
% Dr. Jeremy Kepner
% MIT Lincoln Laboratoy
% kepner@ll.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize MPI.
MPI_Init;

% Create communicator.
comm = MPI_COMM_WORLD;

% Modify common directory from default for better performance.
% comm = MatMPI_Comm_dir(comm,'/tmp');
% comm = MatMPI_Comm_dir(comm,'/gigabit/node-a');

% Get size and rank.
comm_size = MPI_Comm_size(comm);
my_rank = MPI_Comm_rank(comm);

% Do a synchronized start.
starter_rank = 0;
delay = 30;  % Seconds
synch_start(comm,starter_rank,delay);

% Set image size (use powers of 2).
% n_image_x = 2.^17;
% n_image_x = 2.^12;
n_image_x = 2.^(10+1)*comm_size;
n_image_y = 2.^10;

% Number of points to put in each sub-image.
n_point = 100;

% Set filter size (use powers of 2).
n_filter_x = 2.^5;
n_filter_y = 2.^5;

% Set the number of times to filter.
n_trial = 2;

% Computer number of operations.
total_ops = 2.*n_trial*n_filter_x*n_filter_y*n_image_x*n_image_y;

if(rem(n_image_x,comm_size) ~= 0)
 disp('ERROR: processors need to evenly divide image');
 exit;
end

% Print rank.
disp(['my_rank: ',num2str(my_rank)]);

% Set who is source and who is destination.
left = my_rank - 1;
if (left < 0)
  left = comm_size - 1;
end
right = my_rank + 1;
if (right >= comm_size)
  right = 0;
end


% Create a unique tag id for this message (very important in Matlab MPI!).
tag = 1;

% Create timing matrices.
start_time = zeros(n_trial);
end_time = start_time;

% Get a zero clock.
zero_clock = clock;

% Compute sub_images for each processor.
n_sub_image_x = n_image_x./comm_size;
n_sub_image_y = n_image_y;

% Create starting image and working images..
sub_image0 = rand(n_sub_image_x,n_sub_image_y).^10;
sub_image = sub_image0;
work_image = zeros(n_sub_image_x+n_filter_x,n_sub_image_y+n_filter_y);

% Create kernel.
x_shape = sin(pi.*(0:(n_filter_x-1))./(n_filter_x-1)).^2;
y_shape = sin(pi.*(0:(n_filter_y-1))./(n_filter_y-1)).^2;
kernel = x_shape.' * y_shape;

% Create box indices.
lboxw = [1,n_filter_x/2,1,n_sub_image_y];
cboxw = [n_filter_x/2+1,n_filter_x/2+n_sub_image_x,1,n_sub_image_y];
rboxw = [n_filter_x/2+n_sub_image_x+1,n_sub_image_x+n_filter_x,1,n_sub_image_y];

lboxi = [1,n_filter_x/2,1,n_sub_image_y];
rboxi = [n_sub_image_x-n_filter_x/2+1,n_sub_image_x,1,n_sub_image_y];


% Set start time.
start_time = etime(clock,zero_clock);


% Loop over each trial.
for i_trial = 1:n_trial

  % Copy center sub_image into work_image.
  work_image(cboxw(1):cboxw(2),cboxw(3):cboxw(4)) = sub_image;

  if (comm_size > 1)
    % Create message tag.
    ltag = 2.*i_trial;
    rtag = 2.*i_trial+1;

    % Send left sub-image.
    l_sub_image = sub_image(lboxi(1):lboxi(2),lboxi(3):lboxi(4));
    MPI_Send(  left, ltag, comm, l_sub_image );

    % Receive right padding.
    r_pad = MPI_Recv( right, ltag, comm );
    work_image(rboxw(1):rboxw(2),rboxw(3):rboxw(4)) = r_pad;

    % Send right sub-image.
    r_sub_image = sub_image(rboxi(1):rboxi(2),rboxi(3):rboxi(4));
    MPI_Send( right, rtag, comm, r_sub_image );

    % Receive left padding.
    l_pad = MPI_Recv( left, rtag, comm );
    work_image(lboxw(1):lboxw(2),lboxw(3):lboxw(4)) = l_pad;

  end

  % Compute convolution.
  work_image = conv2(work_image,kernel,'same');
  % Extract sub_image.
  sub_image = work_image(cboxw(1):cboxw(2),cboxw(3):cboxw(4));
end


% Get end time for the this message.
end_time = etime(clock,zero_clock);

% Print the results.
total_time = end_time - start_time

% Print compute performance.
total_ops;
gigaflops = total_ops / total_time / 1.e9;
disp(['GigaFlops: ',num2str(gigaflops)]);


% Write data to a file.
outfile = ['blurimage.',num2str(my_rank),'.mat'];
% save(outfile,'start_time','end_time','total_time','kernel','sub_image','work_image');
% save(outfile,'start_time','end_time','total_time','kernel');

% Finalize Matlab MPI.
disp('SUCCESS');
MPI_Finalize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2002 Massachusetts Institute of Technology
% 
% Permission is herby granted, without payment, to copy, modify, display
% and distribute this software and its documentation, if any, for any
% purpose, provided that the above copyright notices and the following
% three paragraphs appear in all copies of this software.  Use of this
% software constitutes acceptance of these terms and conditions.
%
% IN NO EVENT SHALL MIT BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF
% THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF MIT HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% 
% MIT SPECIFICALLY DISCLAIMS ANY EXPRESS OR IMPLIED WARRANTIES INCLUDING,
% BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.
%
% THIS SOFTWARE IS PROVIDED "AS IS," MIT HAS NO OBLIGATION TO PROVIDE
% MAINTENANCE, SUPPORT, UPDATE, ENHANCEMENTS, OR MODIFICATIONS.

