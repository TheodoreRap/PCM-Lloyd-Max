function [xx,xq_x,yy,xq_y,centers,D] = Lloyd_Max_2D(x,N,min_value,max_value)
xx=zeros(length(x)/2,1); % τα διασπόμενα xx και yy απο το x
yy=zeros(length(x)/2,1);
xq_x=zeros(length(x)/2,1); % τα νέα xq_x και xq_y που προκύπτουν απο 
xq_y=zeros(length(x)/2,1); % τα xx  και yy
layers=2^N; % αριθμός επιπέδων κάθετα και άλλα τόσα οριζόντια
centers={layers,layers}; % κέντα περιοχών κβάντισης layers^2
e=10^-6;

j=1; % διαχωρισμός του x σε xx και yy, φτίαχνω ζεύγη xx(1) μαζί με yy(1)
for i=1:2:length(x)
    xx(j)=x(i);
    yy(j)=x(i+1);
    j=j+1;
end

tempo=(max_value-min_value)/layers; % ίδιος τρόπος αρχικοποίησης κέντρων
% αλλά έχουμε 4 περιπτώσεις, 1η να βρισκόμαστε σε xx>0 και yy>0 δηλαδή 1ο
% τεταρτημόριο, 2η να βρισκόμαστε σε xx> και yy<0 δηλαδή 2ο τεταρτημόριο,
% 3η να βρισκόμαστε σε xx<0 και yy<0 δηλαδή 3ο τεταρτημόριο και τέλος xx<0
% και yy>0 δηλασή 4ο τεταρτημόριο

for i=1:layers % αρχικοποίηση κέντρων βάσει τεταρτημορίων
    for j=1:layers
        centers(i,j)={[min_value+tempo/2+(j-1)*tempo,max_value-tempo/2-(i-1)]};
    end
end

repeats=0;
D=[0 1]; % αρχικοποίηση μέσης παραμόρφωσης
d=2;

while abs(D(d)-D(d-1))>=e 
    repeats=repeats+1; % επαναλήψεις μέχρι να έχουμε σύγκλιση
    sum_deformation=0;
    sum_x=zeros(layers,layers); % άθροισμα όλων των xx και yy που περιέχονται  
    sum_y=zeros(layers,layers);% στο κέντρο i,j και το αντίστοιχο πλήθος
    count_xy=zeros(layers,layers); 
    
    for k=1:length(xx)
        min_k=1000;
        for i=1:length(centers)
            for j=1:length(centers)
                m=sqrt((xx(k)-centers{i,j}(1))^2)+((yy(k)-centers{i,j}(2))^2);
                if min_k>m
                    min_k=m;
                    min_x=i;
                    min_y=j;
                end
            end
        end
        xq_x(k)=centers{min_x,min_y}(1);
        xq_y(k)=centers{min_x,min_y}(2);

        sum_deformation=sum_deformation+min_k;
        sum_x(min_x,min_y)=sum_x(min_x,min_y)+xx(k);
        sum_y(min_x,min_y)=sum_y(min_x,min_y)+yy(k);
        count_xy(min_x,min_y)=count_xy(min_x,min_y)+1;
    end
    
    avg_deformation=sum_deformation/(length(x)/2); % υπολογισμός μέσης 
    D=[D avg_deformation]; % παραμόρφωσης για κάθε επανάληψη
    d=d+1;
    
    for i=1:layers
        for j=1:layers % υπολογισμός των νέων κέντρων βάσει της μέσης τιμής
            if count_xy(i,j)~=0 % των x(i) που βρίσκονται μέσα στα αντίσοιχα όρια
                centers(i,j)={[sum_x(i,j)/count_xy(i,j),sum_y(i,j)/count_xy(i,j)]};
            end
        end
    end
end

D(1)=[]; % διαγραφή των πρώτων δυο τιμών του D που είχαμε αρχικοποιήσει
D(1)=[]; % στο 0 και 1

end
