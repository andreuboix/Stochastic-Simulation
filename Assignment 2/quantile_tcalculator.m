i = 1;
studenttdistribution_quantiles = [];
while i <= 30
    studenttdistribution_quantiles = [studenttdistribution_quantiles, tinv(0.975, i-1)];
    i = i+1;
end

display(studenttdistribution_quantiles)
