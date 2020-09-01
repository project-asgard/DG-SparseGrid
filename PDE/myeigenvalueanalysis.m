DT = zeros(5,4);
largestEig = zeros(5,4);
conditionNum = zeros(5,4); %condition number of I-dt*A, with dt = 1/2^level
DOF = zeros(5,4); %degrees of freedom


for i=4:8 %level between 4 and 8
    dt = (1/2)^i; %choose dt same as dx. Chosen for calculating cond
    for j=2:5 %degree between 2 and 5

        %asgard(advection1,'lev',i,'deg',j,'implicit',true,'num_steps',1,'dt',1e-5)
        %asgard(diffusion1,'lev',i,'deg',j,'implicit',true,'num_steps',1)
        
        %asgard(fokkerplanck1_5p1a_noLHS,'lev',i,'deg',j,'implicit',true,'num_steps',1)
        asgard(fokkerplanck1_4p1a,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)
        
        filename1 = sprintf('matrix_iter%03d.mat', 0);
        A = load(filename1,'A');
        A = A.A;
        A = full(A);
        
        lambda = eig(A);
        lambda = lambda(abs(lambda) > 1e-1);
        %plot(lambda,'ro','LineWidth',2)
       
        abslambda = lambda.*conj(lambda); %(norm of lambda)^2
        
        dtEstimate = (-lambda-conj(lambda))./abslambda; %-2Re(lambda)/|lambda|^2
        dtEstimate = dtEstimate(dtEstimate > 0);
        dtEstimate = min(dtEstimate);
        
        DT(i-3,j-1) = dtEstimate; %degree, level
        
        largestEig(i-3,j-1) = max(abs(lambda)); %largest eignevalue in magnitude
        
        [N,~] = size(A);
        I = eye(N);
        conditionNum(i-3,j-1) = cond(I-dt*A);
        DOF(i-3,j-1) = (j+1)*2^i;
    end
end

%%%%%%%%%%%%%%%% plot level vs log of degree %%%%%%%%%%%%%%%%%%%%
level = 4:8;

deg2 = DT(:,1); %first column of DT is degree 2
deg3 = DT(:,2);
deg4 = DT(:,3);
deg5 = DT(:,4);

TD2 = log(deg2); %deg2 => slope = -1.0212
TD3 = log(deg3);
TD4 = log(deg4); %deg4 => slope = -1.2009
TD5 = log(deg5); %deg5 => slope = -1.2381
LEV = [1 4;1 5;1 6;1 7;1 8];

bD2 = LEV\TD2;
bD3 = LEV\TD3;
bD4 = LEV\TD4;
bD5 = LEV\TD5;

DEG2 = exp(LEV*bD2);
DEG3 = exp(LEV*bD3);
DEG4 = exp(LEV*bD4);
DEG5 = exp(LEV*bD5);


dx = zeros(5,1);
for i=4:8
    dx(i-3) = pi/(2^i);
end

%plot with regression
semilogy(level,deg2, level,deg3, level,deg4, level,deg5, level,DEG2, level,DEG3, level,DEG4, level,DEG5, 'LineWidth',2)

lgd = legend('Deg2','Deg3','Deg4','Deg5','Deg2 reg, slope=-1.02','Deg3 reg, slope=-1.13','Deg4 reg, slope=-1.20','Deg5 reg, slope=-1.24');
lgd.NumColumns = 2;
title('Level vs the log of max dt for forward Euler method')

%plot without regression
semilogy(level,deg2, level,deg3, level,deg4, level,deg5, 'LineWidth',2)

legend('Deg2','Deg3','Deg4','Deg5')
title('Level vs the log of max dt for forward Euler method')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%% lolog plot deg vs dt %%%%%%%%%%%%%%%%
degree = 2:5;

lev4 = DT(1,:); %first row of DT is level 4
lev5 = DT(2,:);
lev6 = DT(3,:);
lev7 = DT(4,:);
lev8 = DT(5,:);

TL4 = log(lev4'); %lev4 => slope = -2.0797
TL5 = log(lev5'); %lev5 => slope = -2.3345
TL6 = log(lev6');
TL7 = log(lev7'); %lev7 => slope = -2.1822
TL8 = log(lev8');
DEG = [1 log(2); 1 log(3); 1 log(4); 1 log(5)];

bL4 = DEG\TL4;
bL5 = DEG\TL5;
bL6 = DEG\TL6;
bL7 = DEG\TL7;
bL8 = DEG\TL8;

LEV4 = exp(DEG*bL4);
LEV5 = exp(DEG*bL5);
LEV6 = exp(DEG*bL6);
LEV7 = exp(DEG*bL7);
LEV8 = exp(DEG*bL8);


%plot with regression
loglog(degree,lev4, degree,lev5, degree,lev6, degree,lev7, degree,lev8, degree,LEV4, degree,LEV5, degree,LEV6, degree,LEV7, degree,LEV8, 'LineWidth',2)
lgd = legend('Lev4','Lev5','Lev6','Lev7','Lev8','Lev4 reg, slope=-2.07','Lev5 reg, slope=-2.32','Lev6 reg, slope=-2.56','Lev7 reg, slope=-2.82','Lev8 reg, slope=-3.04');
lgd.NumColumns = 2;
title('log of degree vs log of max dt for forward Euler method')

%plot without regression
loglog(degree,lev4, degree,lev5, degree,lev6, degree,lev7, degree,lev8, 'LineWidth',2)
legend('Lev4','Lev5','Lev6','Lev7','Lev8')
title('log of degree vs log of max dt for forward Euler method')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating the CFL
CFL = zeros(5,4);

CFL(1,:) = lev4*2^8; %CFL for level 4
CFL(2,:) = lev5*2^10;
CFL(3,:) = lev6*2^12;
CFL(4,:) = lev7*2^14;
CFL(5,:) = lev8*2^16;
% B = 0;
% Q = 0;
% for d=2:5
%     for n=4:8
%         B = B+DT(n-3,d-1)*d^-2.56*exp(-1.15*n);
%         Q = Q+(d^-2.56*exp(-1.15*n))^2;
%     end
% end
% B = B/Q;
% 
% N = zeros(4,5);
% for d=2:5
%     for n=4:8
%         N(d-1,n-3) = 19.2*d^(-2.56)*exp(-1.15*n);
%     end
% end
% 
% norm(N-DT');