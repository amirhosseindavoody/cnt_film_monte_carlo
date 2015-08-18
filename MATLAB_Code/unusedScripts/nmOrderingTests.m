close all;

i=1;
for n=1:100
    for m=0:n
        out(i) = n*100+m;
        i = i + 1;
    end
end

count = zeros(1,max(out));

for i=1:length(out)
   
    count(out(i)) = count(out(i)) + 1;
    
end

figure;
bar(count);