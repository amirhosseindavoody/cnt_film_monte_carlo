% Simply looking at the plots of the random numbers vs the possible trs
% that can come out of it.


% plot of possible tr's 

close all;

gamma = 2;
r = linspace(0,1,1000);

tr = -(1/gamma)*log(r);


plot(r, tr);


exp(1/gamma*gamma*-1)
