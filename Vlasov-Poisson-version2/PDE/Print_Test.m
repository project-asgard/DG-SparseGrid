function[a] = re_print_functions()
  ne = 1e19;
  Te = 10;
  Zeff = 1.4;
  v = 1;
  ve = re_vThermal(Te)
  x = re_x(v,ve)
  Psi = re_Psi(x)
  coulombLog = re_coulombLog(ne,Te)
  Gamma = re_Gamma(ne,coulombLog)
  delta = re_delta(ve)
  vee = re_vee(ne,Te,ve,coulombLog)
  Ec = re_Ec(ne,coulombLog)
  lowerGamma = re_lowerGamma(v)
  rho = re_rho(lowerGamma,v)
  C_A = re_C_A(Gamma,v,Psi)
  C_B = re_C_B(Zeff,Gamma,delta,v,x,Psi)
  C_F = re_C_F(Te,Gamma,Psi)
  a = 0;
endfunction