
function ran = RandomAccess_rand(ran);
% Create random 64 bit numbers for
% RandomAccess benchmark.
 
  % Define constants used in random number generator.
  POLY = uint64(7);

  % Get blocksize.
  BLOCKSIZE = length(ran);

  % Vectorized approach.
  % Get elements in ran that set 64th bit set.
  i_bit64 =find(bitget(ran,64));

  % Do shift.              
  ran = bitshift(ran,1);
      
  % Apply XOR.
  if not(isempty(i_bit64))
    ran(i_bit64) =  bitxor(ran(i_bit64),POLY);
  end
