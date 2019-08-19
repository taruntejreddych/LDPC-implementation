clear
clc
snr=[1.2 1.4 1.6 1.8 2];
for value=1:length(snr)
iter = 1 ;
var(value) = 10^(-(snr(value)/10));
data=load('h.mat');
M=sparse(data.H);
parity=comm.LDPCEncoder(M);
H_matrix = full(data.H);

Img = imread('bitmap.bmp')/255;
inp=reshape(Img,[2000,1]);
inp=double(inp);
%inp = logical(Img.');

%H=mod(rref(data.H),2);
% inp=logical(randi([0 1], 2000, 1));
% enc=comm.LDPCEncoder(inp);
encode=step(parity,inp);
encode=transpose(encode);

modulated_message=(-1).^encode;

noise = sqrt(var(value))*randn(1,length(modulated_message));
recieved_message = modulated_message + noise;
LLR = recieved_message*2/var(value);

kappa = zeros(length(H_matrix(:,1)),length(H_matrix(1,:)));
lambda = zeros(length(H_matrix(1,:)),length(H_matrix(:,1)));

for count = 1:length(H_matrix(1,:))
    ones_index = find(H_matrix(:,count));
    lambda(count,ones_index) = LLR(count);
end

%for iter = 1:num_iterations
while iter<10
    
    for i=1:length(H_matrix(:,1))
        ones=find(H_matrix(i,:));
        product=1;
        for k=1:length(ones)
            product=product*tanh(lambda(ones(k),i)/2);
        end
        product_vect=product./tanh(lambda(ones,i)/2);
        kappa(i,ones)=2*atanh(product_vect);
    end
    
    for count = 1:length(H_matrix(1,:))
        ones_index = find(H_matrix(:,count));
        lambda(count,ones_index) = LLR(count);
    end
    llr_f = zeros(length(H_matrix(1,:)));
    for j = 1:length(H_matrix(1,:))
        ones_index = find(H_matrix(:,j));
        for k = 1:length(ones_index)
            llr_f(j) = llr_f(j) + kappa(ones_index(k),j);
            for i=1:length(ones_index)
                if (k~=i)
                    lambda(j,ones_index(k)) = lambda(j,ones_index(k)) + kappa(ones_index(i),j);
                end
            end
           % lambda(j,ones_index(k)) = lambda(j,ones_index(k)) + LLR(j);
        end
        llr_f(j) = llr_f(j) + LLR(j);
    end
    

    c = llr_f;
    demodulated_message = zeros(1,length(encode));
    for u = 1:length(encode)
        if c(u) > 0
        demodulated_message(u) = 1;
        stop_msg(u) = 1;
        end
        if c(u)<0
            demodulated_message(u) = -1;
            stop_msg(u) = 0;
        end
    end
    
    
   % LLR = llr_f;
%     stop_condition = ~any(mod(H_matrix*transpose(demodulated_message),2))
     %stop = logical(double(H_matrix)*double(transpose(stop_msg)));
     stop = binary_multiplication(full(H_matrix),transpose(stop_msg));
  
    if( isequal(stop,zeros(2000,1)))
        iter
        break
    end
    if (mod(iter,2)==0)
    %iter
    stop = length(find(stop))
    bits_flipped(value) = length(find(modulated_message-demodulated_message))
    end
   iter = iter+1;
end
bits_flipped(value) = length(find(modulated_message-demodulated_message))
a = comm.LDPCDecoder(M);
q = step(a,transpose(stop_msg));
%[n,r] = biterr(encode,stop_msg)
subplot(2,length(snr),value);
imshow(reshape(q,[40,50]));
title(['Using LDPC and snr =',num2str(snr(value))]);

end
% scatter(snr,bits_flipped/4000);
% line(snr,bits_flipped/4000);
% set(gca,'xscale','log')
% semilogy(snr,bits_flipped/4000,'o');
% line(snr,bits_flipped/4000);
% grid on;
%------------------------------------------------------------
for value=1:length(snr)
data=load('h.mat');
M=sparse(data.H);
parity=comm.LDPCEncoder(M);
H_matrix = data.H;

Img = imread('bitmap.bmp')/255;
inp=reshape(Img,[2000,1]);
inp=double(inp);

encode=step(parity,inp);
encode=transpose(encode);

modulated_message=(-1).^encode;

noise = sqrt(var(value))*randn(1,length(modulated_message));
recieved_message = modulated_message + noise;

c = recieved_message;
demodulated_message = zeros(1,length(encode));
for u = 1:length(encode)
    if c(u) > 0
    demodulated_message(u) = 1;
    stop_msg(u) = 1;
    end
    if c(u)<0
        demodulated_message(u) = -1;
        stop_msg(u) = 0;
    end
end
bits_flipped = length(find(modulated_message-demodulated_message))
dem = comm.LDPCDecoder(M);
q = step(dem,transpose(stop_msg));
%[n,r] = biterr(encode,stop_msg)
subplot(2,length(snr),length(snr)+value);
imshow(reshape(q,[40,50]));
title(['Without LDPC and snr=',num2str(snr(value))]);
end