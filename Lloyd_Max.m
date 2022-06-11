function [xq,centers,D] = Lloyd_Max(x,N,min_value,max_value)
xx=zeros(length(x),1);
xq=zeros(length(x),1);
layers=2^N; % αριθμός επιπέδων
centers=zeros(layers,1); % κέντα περιοχών κβάντισης
e=10^-6; 
limits=zeros(layers+1,1); % όρια

tempo=(max_value-min_value)/layers; 

for i=1:layers % επιλογή αρχικών κέντρων στην μέση κάθε ορίου
    centers(i)=min_value+tempo/2+(i-1)*tempo;
end

repeats=0;
D=[0 1]; % αρχικοποίηση μέσης παραμόρφωσης
d=2;

while abs(D(d)-D(d-1))>=e 
    repeats=repeats+1; % επαναλήψεις μέχρι να έχουμε σύγκλιση
    sum_deformation=0;
    
    limits(1)=min_value; % η ελάχιστη τιμή είναι το πρώτο όριο
    for i=2:layers % υπολογισμός λοιπών ορίων κρατόντας σταθερά τα άκρα
        limits(i)=(centers(i-1)+centers(i))/2; % των ορίων
    end
    limits(i+1)=max_value; % η μέγιστη τιμή είναι το τελευταίο όριο

    sum_x=zeros(layers,1); % άθροισμα όλων των x(i) που περιέχονται στο 
    count_x=zeros(layers,1); % κέντρο i και το αντίστοιχο πλήθος

    for i=1:length(x)
        for j=1:(length(limits)-1)
            if x(i)>limits(j) && x(i)<limits(j+1) % έλεγχος για την ακριβής
                xx(i)=centers(j); % περιοχή του x(i)
                xq(i)=layers+1-j; % τιμές 1:2^Ν με 1 στο μεγαλύτερο layer
                sum_deformation=sum_deformation+abs(centers(j)-x(i));
                sum_x(j)=sum_x(j)+x(i);
                count_x(j)=count_x(j)+1;
            end
            
        end
        if x(i)==min_value % περίπτωση που το x(i) έχει την ελάχιστη
            xx(i)=centers(1);
            xq(i)=layers;
            sum_deformation=sum_deformation+abs(centers(1)-x(i));
            sum_x(1)=sum_x(1)+x(i);
            count_x(1)=count_x(1)+1;
        end
        if x(i)==max_value % ή την μέγιστη τιμή
            xx(i)=centers(layers);
            xq(i)=1;
            sum_deformation=sum_deformation+abs(centers(layers)-x(i));
            sum_x(layers)=sum_x(layers)+x(i);
            count_x(layers)=count_x(layers)+1;
        end
    end

    avg_deformation=sum_deformation/length(x); % υπολογισμός μέσης 
    D=[D avg_deformation]; % παραμόρφωσης για κάθε επανάληψη i του
    d=d+1; % αλγορίθμου

    for j=1:layers % υπολογισμός των νέων κέντρων βάσει της μέσης τιμής των
        if count_x(j)~=0 % x(i) που βρίσκονται μέσα στα αντίσοιχα όρια
            centers(j)=sum_x(j)/count_x(j);
        end
    end
end

D(1)=[]; % διαγραφή των πρώτων δυο τιμών του D που είχαμε αρχικοποιήσει
D(1)=[]; % στο 0 και 1

end
