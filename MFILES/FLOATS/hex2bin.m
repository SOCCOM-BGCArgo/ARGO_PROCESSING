function s = hex2bin(x,n)
%HEX2BIN Convert a hexadecimal string into a binary string
% str = hex2dec(x) produces a binary representation of x as a string.
%   The value of x must be smaller than hexadecimal 0x10000000000000.
% str = hex2dec(x, n) produces a binary representation with at least n bits.
%
% The output of hex2bin is independent of the endian settings
% of the computer you are using.
%
% See also DEC2BIN, BIN2DEC, DEC2HEX, HEX2DEC

if nargin<2
    s = dec2bin(hex2dec(x));
else
    s = dec2bin(hex2dec(x),n);
end
end