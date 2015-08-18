%% Builds test rate tables

clear all;
close all;

%Chiralities
chirList = Chirality(5,1);

chirList(1).n = 7;
chirList(1).m = 5;

chirList(2).n = 7;
chirList(2).m = 6;

chirList(3).n = 8;
chirList(3).m = 6;

chirList(4).n = 8;
chirList(4).m = 7;

chirList(5).n = 9;
chirList(5).m = 7;

size = 10;

rtTable = zeros(size,4*size,length(chirList)^2); 

chirCount = 1;
for i=1:length(chirList)
    for j=1:length(chirList)
       for k=1:size
        for l=1:4*size
            rtTable(k,l,chirCount) = k*chirList(i).n*3+l*chirList(i).m^(3/2) + chirList(j).m;
        end
       end
      chirCount = chirCount + 1;
    end
end

chirCount = 1;
for i=1:length(chirList)
   for j=1:length(chirList)
       filename = ['output\' num2str(chirList(i).n) ',' num2str(chirList(i).m) '_' ...
           num2str(chirList(j).n) ',' num2str(chirList(j).m) '.bin'];
       file = fopen(filename,'w');
       fwrite(file,rtTable(:,:,chirCount),'double');
       fclose(file);
       chirCount = chirCount + 1;
   end
end


file = fopen('output\7,5_7,5.bin');
x = fread(file,[size,4*size],'double');



