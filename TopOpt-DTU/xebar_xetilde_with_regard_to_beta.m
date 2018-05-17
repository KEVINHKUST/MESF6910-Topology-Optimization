n = 1;
b = 1;
i = 1;
while(b <= 512)
    figure(n)
    xetilde = 0:0.001:1;
    xebar = 1 - exp(-b*xetilde) + xetilde * exp(-b);
    plot(xetilde,xebar,'color',[rand,rand,rand])
    str{i}=['¦Â =',num2str(b)];
    hold on
    b = b * 2 ;
    i = i + 1 ;
end
legend(str,'location','best');
ylabel('xe__bar')
xlabel('xe__tilde')


