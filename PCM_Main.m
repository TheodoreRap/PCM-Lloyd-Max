%------------------------source Α------------------------------------------
x=double(randn(10000,1));
min_value_x=x(1);
max_value_x=x(1);
for i=2:length(x)
    if min_value_x>x(i)
        min_value_x=x(i);
    end
    if max_value_x<x(i)
        max_value_x=x(i);
    end
end

%------------------------source B------------------------------------------
b=1;
a=[1 1/2 1/3 1/4 1/5 1/6];
y=filter(b,a,x); % τυχαία διαδικασία AR 5ης τάξης
min_value_y=y(1);
max_value_y=y(1);
for i=2:length(y)
    if min_value_y>y(i)
        min_value_y=y(i);
    end
    if max_value_y<y(i)
        max_value_y=y(i);
    end
end

%---------------------------N=2,3,4----------------------------------------
% κάθε φορά αλλάζουμε το Ν=2,3,4 και για την πηγή Α: x, min_value_x,
% max_value_x και την πηγή Β: y, min_value_y, max_value_y
N=2;
[xq,centers,D] = Lloyd_Max(x,N,min_value_x,max_value_x);
bit_xq=de2bi(xq-1,'left-msb'); % δυαδική αναπαράσταση N bits

[xx,xq_x,yy,xq_y,centers_2,D_2] = Lloyd_Max_2D(x,N,min_value_x,max_value_x);

%-----------------------------SQNR-----------------------------------------
SQNR=zeros(length(D),1); % υπολογισμός SQNR βαθμωτού κβαντιστή
for i=1:length(D)
    SQNR(i)=10*log10(mean(x.^2/D(i))); % SQNR σε dB
end

SQNR_2=zeros(length(D_2),1); % υπολογισμός SQNR διανυσματικού κβαντιστή
for i=1:length(D_2)
    SQNR_2(i)=10*log10(mean(((xx.^2+yy.^2)/2)/D_2(i))); % SQNR σε dB
end

%------------------------------MSE-----------------------------------------
centers=flip(centers); % ώστε να είναι σωστή η αντιστοιχία 1,2,..,2^Ν
z=centers(xq);
MSE=0; % μέσο τετραγωνικό σφάλμα (MSE) βαθμωτού
for i=1:length(x)
    MSE=MSE+ (z(i)-x(i))^2;
end
MSE=MSE/length(x);

MSE_2=0; % μέσο τετραγωνικό σφάλμα (MSE_2) διανυσματικού
for i=1:length(xx)
    MSE_2=MSE_2+ (sqrt((xx(i)-xq_x(i))^2)+((yy(i)-xq_y(i))^2))^2;
end
MSE_2=MSE_2/length(xx);


%--------------------------plots-scatters----------------------------------

% plot που περιέχει τα SQNR - REAPEATS και D - REPEATS
% plot(SQNR)
% plot(SQNR_2)
% plot(D)
% plot(D_2)

% plot που περιέχει την είσοδο x ή y και το xq του βαθμωτού κβαντιστή
% hold on
% plot(x)
% plot(z)
% hold off

% scatter που περιέχει την είσοδο διάνυσμα xx και yy καθώς και το xq_x και 
% το xq_y του διανυσματικού κβαντιστή
% hold on
% scatter(xx,yy)
% scatter(xq_x,xq_y)
% hold off
