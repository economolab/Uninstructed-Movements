weight = [0.7961, 0.5705, 0.2929, 0.1504, 0.0772, 0.0204];
times = [0.3,0.6,1.2,1.8,2.4,3.6];
bar(times,weight,0.4,"grouped",'black')
x = linspace(0.3,3.6);
y = exppdf(x,0.9);
hold on; 
plot(x,y,'LineWidth',2,'Color','green')
xlabel('Delay length (s)')
ylabel('PDF (prob of delay = x)')