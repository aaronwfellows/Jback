function r2 = rsq(A,B)
corrcoefficient=corrcoef(A,B);
if size(corrcoefficient)~= [2 2]
   r2 = NaN;
else
   r2=(corrcoefficient(1,2))^2;
end
return
