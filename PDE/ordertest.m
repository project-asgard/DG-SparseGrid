%With source
deg1Error = [4.1746e-02 2.1496e-02 3.3106e-02 4.4622e-02];
deg2Error = [3.7094e-02 1.0947e-02  2.9150e-03 7.1458e-04  1.7441e-04 4.3154e-05];
deg3Error = [9.4304e-04 1.1529e-04 2.1238e-05 3.6916e-06 5.8426e-07 8.6636e-08];

%without source

order = zeros(1,5);
for i=1:5
    order(i) = (log(deg3Error(i)/deg3Error(i+1)))/log(2);
end
order