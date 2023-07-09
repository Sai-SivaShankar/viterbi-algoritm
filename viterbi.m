clc;
clear;
close all;
%% Generating a random bit string of size N=256 
N = 256;
databits = randi([0,1], 1, N);
D = [0 0 databits 0 0];
%% Encoding
g1 = [1 1 0];
g2 = [1 1 1];
g3 = [1 0 1];
C = zeros(1,length(D)*3);
reg = zeros(1,3);
for i = 1:length(D)-2
    reg(1:3) = D(i:i+2);
    C(3*i-2) = xor(reg(3),reg(1));
    C(3*i-1) = xor(reg(1),xor(reg(2),reg(3)));
    C(3*i) = xor(reg(1),reg(2));
end
C = C(1:3*(N+2));
%% trellis

%A trellis should specify the nextstates and output for a given input and 
%current state.
numInputs = 2;
numOutputs = 8;
numStates = 4;

nextState = zeros(numStates, numInputs);%It is a numStates x numInputs size matrix. 
outputs = zeros(numStates, numInputs);%It is a numStates x numInputs size matrix.

for k=1:numStates%current_state
    nextState(k,1) = floor((k-1)/2);
    nextState(k,2) = 2 + floor((k-1)/2);
    
    outputs(k,1) = pow2(2)*xor(xor(g1(3)*0, g1(2)*floor((k-1)/2)),g1(1)*(mod(k-1,2))) + pow2(1)*xor(xor(g2(3)*0, g2(2)*floor((k-1)/2)),g2(1)*(mod(k-1,2))) + pow2(0)*xor(xor(g3(3)*0, g3(2)*floor((k-1)/2)),g3(1)*(mod(k-1,2)));
    outputs(k,2) = pow2(2)*xor(xor(g1(3)*1, g1(2)*floor((k-1)/2)),g1(1)*(mod(k-1,2))) + pow2(1)*xor(xor(g2(3)*1, g2(2)*floor((k-1)/2)),g2(1)*(mod(k-1,2))) + pow2(0)*xor(xor(g3(3)*1, g3(2)*floor((k-1)/2)),g3(1)*(mod(k-1,2)));
end

table(nextState,outputs);

%% Recieved Data
p = 0:0.005:0.5;%probability of error
rsig = zeros(length(p),length(C));
%% Decoding : Vetirbi algorithm
numerrbits = zeros(1,length(p));
biterr_rate = zeros(size(numerrbits));
decData = zeros(length(p),length(databits));
biterr_rateavg = zeros(size(biterr_rate));
for j = 1:10
for i= 1:length(p)
rsig(i,:) = bsc(C, p(i));
decData(i,:) = viterbi(rsig(i,:), numStates, N, nextState, outputs);
[numerrbits(i), biterr_rate(i)] = biterr(decData(i,:), databits);
end
biterr_rateavg = biterr_rate./10 + biterr_rateavg;
end
%% plotting bit error rate vs crossover probability
figure;
plot(p,biterr_rateavg);
xlabel('crossover probability');
ylabel('biterror rate');
%% part B
p1 = 0:0.01:0.5;
p0 = 0:0.01:0.5;

rsig1 = zeros(length(p1),length(C));
biterr_rate1 = zeros(length(p0),length(p1));
for i = 1:length(p0)
    for j = 1: length(p1)
        avg = 0;
        for runs = 1 :4
            rsig1 = basc(C,p0(i),p1(j));
            decData1 = viterbi(rsig1, numStates, N, nextState, outputs);
            [~, u] = biterr(decData1, databits);
            avg = avg + u;
        end
        biterr_rate1(i,j) = avg/10;
    end 
end

%% functions
function op_data = basc(inp, p0, p1)
    op_data = zeros(1,length(inp));
    for k=1:length(inp)
        if inp(k) == 1
            op_data(k) = bsc(inp(k), p1);
        else
            op_data(k) = bsc(inp(k), p0);
        end
    end
end

function decoded = viterbi(rsig, numStates, N, nextState, outputs)
    
    dparr = ones(numStates, N+1)*inf;
    dparr(1,1) = 0;
    %stores the previous element.
    parent = ones(numStates, N+1)*inf;
    parent(1,1) = 0;
    
    for k=1:N
        for j=1:numStates
            if dparr(j,k) ~= inf
                
                dist1 = dparr(j, k) + hammDist(bi2de(rsig(3*k-2:3*k)), outputs(j,1));
                dist2 = dparr(j, k) + hammDist(bi2de(rsig(3*k-2:3*k)), outputs(j,2));
                
                if dist1 < dparr(nextState(j,1)+1, k+1)
                    dparr(nextState(j,1)+1, k+1) = dist1;
                    parent(nextState(j,1)+1, k+1) = j;
                end
                if dist2 < dparr(nextState(j,2)+1, k+1)
                    dparr(nextState(j,2)+1, k+1) = dist2;
                    parent(nextState(j,2)+1, k+1) = j;
                end
                
            end
        end
    end
    ind = N+1;
    [~,state]  = min(dparr(:, N+1));
    decoded = zeros(1,length(rsig)/3);
    while ind > 1
        decoded(ind) = floor((state-1)/2);
        state = parent(state,ind);
        ind = ind-1;
    end
    decoded = decoded(2:end-1);
end

function hamm_dist = hammDist(num1, num2)
    temp2 = bitxor(num1, num2);
    temp2 = de2bi(temp2);
    hamm_dist = sum(temp2 == 1);
end





