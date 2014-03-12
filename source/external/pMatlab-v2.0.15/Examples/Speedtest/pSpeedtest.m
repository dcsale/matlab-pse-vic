%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script times SendMsg/RecvMsg for
% a variety of message sizes.
% To run in parallel 
% at the Matlab prompt type 
%   eval(pRUN('pSpeedtest',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nmessage = 20;  % Number of message sizes.
Ntrial = 4;   % Trials per message size.

if(Np < 2)
 disp('ERROR: too few processors (need at least 2)');
 exit;
end

source = mod(Pid - 1,Np);  % Set source.
dest = mod(Pid + 1,Np);   % Set destination.

tag = 8;  % Initial messge tag.
TotalTime = zeros(Ntrial,Nmessage);  % Timing matrix.

% Compute message sizes.
p = 1:Nmessage;
messageSize = 2.^p;
ByteSize = 8.*messageSize;

for i_message = 1:Nmessage  % Loop over each message size.

  % Create message.
  sendData = zeros(1,messageSize(i_message)) + Pid;

  for i_trial = 1:Ntrial

    tic; % Start clock.
    SendMsg(dest,tag,sendData);  % Send data.
    recvData = RecvMsg(source,tag);  % Receive data.
    TotalTime(i_trial,i_message) = toc;  % Stop clock.

    % Check data.
    if(any(recvData ~= source))
      disp('WARNING: incorrect data sent.');
    end

    tag = mod(tag + 1,32)+1; % Increment message tag.

  end
end

% Compute bandwidth.
MessageBytes = repmat(ByteSize,Ntrial,1);
Bandwidth = 2.*MessageBytes./TotalTime;

% Write data to a file.
outfile = ['speedtest.',num2str(Pid),'.mat'];
save(outfile,'ByteSize','TotalTime','MessageBytes','Bandwidth');
