function [ S ] = powerlaw(p,Extra)

  global idx;
  format compact;

  Lc = 0.4825;   % borehole radius (m)
  mu = 1.75E-03; % brine viscosity (Pa*sec)
  cf = 3.1E-10; % brine compressibility (1/Pa)
  cm = 10.0^p(5);  % rock compressibility (1/Pa)

  lamda = -999.9; 
  omega = -999.9;
  
  k0 = 10.0.^p(1);
  n0 = 10.0.^p(2);
  eta = p(3);
  tau = p(4);
  m = p(6);
  c = (n0*cf + cm);

  Tc = n0*c*Lc*Lc*mu/k0;  % characteristic time
  Pc = 3.0E+6; % difference between far-field and near-field pressure 
  
  kappa = eta * tau;

  all_time = Extra.time;
  nt_vec = Extra.ndata;
  cmd = './drive-powerlaw.sh';

  S = [];

  % loop over individual observation wells
  for i = 1:4

    if i == 1
        idx0 = 1;
        idx1 = nt_vec(1);
    else
        idx0 = 1 + sum(nt_vec(1:i - 1));
        idx1 = sum(nt_vec(1:i));
    end

    %[idx0,idx1]

    tD = all_time(idx0:idx1)/Tc; % dimensionless time
    save('tD.in','tD','-ascii'); 

    pout = [eta; kappa; m; Extra.rvec(i)]; 
    save('parameters.in', 'pout', '-ascii');

    unix(cmd);
   
    X = load('powerlaw.out');
    S = [S; X(:,2)*Pc]; % p (dimensionless) -> Pa

  end
   
  if mod(idx,20) == 0
    fprintf('%5i %.2f %.2f %.2f %.2f %.2f %.2f\n', ...
              idx,p(1),p(2),p(3),p(4),p(5),p(6));
  end 
   
  idx = idx + 1;
  

  
