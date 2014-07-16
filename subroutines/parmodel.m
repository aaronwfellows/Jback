function nee=parmodel(abc,par)
  
nee = abc(1) + (abc(2)*par)./(par+abc(3));  
  
  