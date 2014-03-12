
function ran = RandomAccessStarts(n);
% Utility routine to start random number generator at Nth step

  % Define constants used in random number generator.
  % PERIOD = uint64(1317624576693539401);
  PERIOD = bitset(bitset(bitset(uint64(1317624576693539401),7),4),1);
  POLY = uint64(7);

  % Get the number of locations that we are starting.
  BLOCKSIZE = length(n);

  % Allocate variables.
  m2 = zeros(1,64,'uint64');
  temp = uint64(1);

  ran = zeros(1,BLOCKSIZE,'uint64');
  ran(:) = uint64(2);

  % UNLIKELY WE WILL ENCOUNTER THESE LIMITS SO IGNORE
  % % Modulate
  % while (n < 0) n = n + PERIOD;
  % while (n > PERIOD) n -= PERIOD;

  % Build up m2 table.
  for i=0:63

    % Initialize m2.
    m2(i+1) = temp;

    % temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
    % temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);

    % Equivalent to above 2 lines of C code.
    for j=1:2
      % Check if bit 64 is set.
      if bitget(temp,64)
        temp = bitshift(temp,1);
        temp = bitxor(temp,POLY);
      else
        temp = bitshift(temp,1);
      end
    end

  end

  % Allocate indices.
  ii = zeros(1,BLOCKSIZE);

  % Initialize index table ii for each starting point.
  % Use "vectorized" method.
  for j=62:-1:0
    % Shift n by i, and check first bit,
    % then find all instance where this statement is true.
    i_hit =find( bitand(bitshift(uint64(n),-j),uint64(1)) );
    if not(isempty(i_hit))
      % Set index to first i_hit.
      ii(i_hit) = max(ii(i_hit),j);
    end
  end

  % Compute ran for these  starting points.
  % To do this in a "vectorized" manner we
  % replace if statements with finds.

  temp = zeros(1,BLOCKSIZE,'uint64');

  % Find subset jj where ii > 0.
  jj = find(ii > 0);
  while ( not(isempty(jj)) )

    % Initalize jj subset temp.
    temp(jj) = uint64(0);


    for j=0:63

      % Left shift ran by j and check first bit,
      % then find sub-subset where this is true.
      jjj = find( bitand(bitshift(ran(jj),-j),uint64(1)) );
      if (not(isempty(jjj)) )

        % XOR the sub-subset jj(jjj) by m2.
        temp(jj(jjj)) = bitxor(temp(jj(jjj)),m2(j+1));

      end
    end 

    % Copy subset to ran and decrement ii.
    ran(jj) = temp(jj);
    ii(jj) = ii(jj) - 1;

    % Left shift n by ii over the set jj
    % then find the sub-subset where the first bit is true.
    jjj = find( bitand(bitshift(uint64(n(jj)),-ii(jj)),uint64(1)) );
    if ( not(isempty(jjj)) )

      % Find the sub-sub-subset where the 64th bit is set.
      jjjj = find( bitget(ran(jj(jjj)),64) );

      % Right shift the sub-subset jj(jjj)  by 1.
      ran(jj(jjj)) = bitshift(ran(jj(jjj)),1);
      if ( not(isempty(jjjj)) )

        % XOR the sub-sub-subset jj(jjj(jjjj)) with POLY.
        ran(jj(jjj(jjjj))) = bitxor(ran(jj(jjj(jjjj))),POLY);

      end

    end

    % See if there and more left to compute.
    jj = find(ii > 0);

  end

  % Set all cases where n = 0 to 1.
  ran( find(n == 0) ) = uint64(1);
