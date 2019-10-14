x=32:16:320;%x轴上的数据，第一个值代表数据开始，第二个值代表间隔，第三个值代表终止
 a=RATE(1:8:145); %a数据y值
 b=32:2:256; %b数据y值
 %plot(x,RS(1:113),'-*b',x,RD(1:113),'-or'); %线性，颜色，标记
 plot(x,a,'-r.'); %线性，颜色，标记
axis([32,320,1,1.5])  %确定x轴与y轴框图大小
set(gca,'XTick',[32:20:320]) %x轴范围1-6，间隔1
set(gca,'YTick',[1:0.1:1.5]) %y轴范围0-700，间隔100
legend('Iteration Ratio Sys/Ring');   %右上角标注
xlabel('Matrix size')  %x轴坐标描述
ylabel('iteration time') %y轴坐标描述