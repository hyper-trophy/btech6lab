data = round(rand(10, 1));

[n m] = size(data);

expanded_data = repmat(data, [1, 8]);

H = commsrc.pn(...
 'GenPoly', [3 2 0], ...
 'InitialStates', [0 0 1], ...
 'CurrentStates', [0 0 1], ...
 'Mask', [0 0 1], ...
 'NumBitsOut', 8);

pn=generate(H); 
pn_repeated = repmat(reshape(pn, [1 numel(pn)]), [n 1]);
xored = xor(expanded_data, pn_repeated);

M = 2;

data_vector = reshape(xored, [numel(xored)/M M]);

decimal_data_vector = bi2de(data_vector);

transmitted = pskmod(decimal_data_vector, 2^M);

SNR = 2;

received = awgn(transmitted, SNR);

demodulated = pskdemod(received, 2^M);

binary_data_vector = de2bi(demodulated);

binary_data = reshape(binary_data_vector, size(expanded_data));

despread_data = xor(binary_data, pn_repeated);

msg_rx = round(mean(despread_data, 2));

BER = mean(abs(msg_rx - data));

fprintf('\nsent data : %s\nreceived data : %s\nBER : %s\n',... 
    num2str(data), num2str(msg_rx), num2str(BER))

