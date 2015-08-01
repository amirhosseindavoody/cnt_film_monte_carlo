%% Testing methods to end simulation

%Could do something like this and if the output value is less than a
%threshold, say .01 for some number of steps, say 5, then the simulation
%will end. 
figure;
excitonCountDeriv = diff(excitonCount);
output = zeros(length(excitonCountDeriv),1);
maxDeriv = 0;
for i=1:length(excitonCountDeriv)
   if(excitonCountDeriv(i) > maxDeriv)
       maxDeriv = excitonCountDeriv(i);
   end
   output(i) = excitonCountDeriv(i)/maxDeriv;
end

plot(sparset(1:end-1),output)

% Running average method
bufSize = 1000;
buf = zeros(1,bufSize);
for i=1:bufSize
    buf(i) = sum(excitonDist(i,:)); 
end

bufSum = sum(buf);
output2 = zeros(1,length(excitonDist)-bufSize);

for i=bufSize+1:length(excitonDist)
    bufSum = bufSum - sum(excitonDist(i-bufSize,:))/bufSize;
    bufSum = bufSum + sum(excitonDist(i,:))/bufSize;
    output2(i-bufSize) = bufSum;
end

plot(output2)