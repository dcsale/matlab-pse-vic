function MatMPI_Sleep
  if (isunix)
    pause(0.01);
  else
    pause(0.1);
  end
  return;
