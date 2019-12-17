k = [];
n = 20;
for i = 1:n
    T = readtable(strcat(num2str(i), '.csv'));
    A = table2array(T);
    x = A(2:1500, 1);
    y = A(2:1500, 2);
    yy = smoothdata(y, 'gaussian'); 
    s = num2str(i);
    f = yy;
    k = [k,f];
    disp(i)
    plot (f)
    hold on
end
M = mean(k,2);
smoothM = smoothdata(M, 'gaussian');
disp(smoothM)
plot(smoothM);
hold off
plot(smoothM);
hold on
p = poly(M);

%sz = size(M);
%disp(sz)
%xn = linspace(0,2,3999);
%xx = xn.';
%p = polyfit(xx,M,100);
%plot (p)
hold off









