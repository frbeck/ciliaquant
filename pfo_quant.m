%Make an array to hold pixel intensity values
k = [];
%Indicate the sample size
n = 2;
%Make a matrix to holld y length values
ylTable = [];
%For each image identify the y length and insert it into the array ylTable
for i = 1:n
    T = readtable(strcat(num2str(i), '.csv'));
    A = table2array(T);
    sz = size(A);
    sz2 = sz(1);
    ylTable = [ylTable, sz2];
end
%Select the smallest y length from ylTable
ymin = min(ylTable);
%Compile the data and normalize against the noise
for i = 1:n
    T = readtable(strcat(num2str(i), '.csv'));
    A = table2array(T);
    x = A(2:ymin, 1);
    y = A(2:ymin, 3);
    y1 = A(1, 3);
    Nr =  y./y1; 
    s = num2str(i);
    f = Nr;
    %insert the normalized data into a new matrix
    k = [k,f];
    disp(i)
    hold on
end
%Average the all of the data arrays
M = mean(k,2);
smoothy = smoothdata(M, 'rloess');
%Smooth the averaged data
smoothM = smoothdata(M, 'gaussian');
%Plot the results
derivM = diff(smoothM);
scatter(x, smoothM)
hold on
plot(M)
%plot(f)
hold off






