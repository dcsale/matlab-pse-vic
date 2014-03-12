function MatMPI_Sleep

% Turn off warning that occurs when running Matlab on Linux with no
%     stdin source attached
%warning('off', 'MATLAB:typeaheadBufferOverflow');
if (isunix)
    pause(0.01);
else
    pause(0.1);
end
return;
