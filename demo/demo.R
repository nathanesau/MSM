# Constant force mortality

print(ax(cfm, n=10, x=60, h=1/12, delta=log(1.05)))
print(Ax(cfm, n=10, x=60, h=1/12, delta=log(1.05)))

tmp1 <- 0.05/(0.05 + log(1.05)) * (1 - exp(-(log(1.05) + 0.05)*10))
print(tmp1)
tmp2 <- (1 - exp(- (log(1.05)+0.05)*10)) / (0.05 + log(1.05))
print(tmp2)

# Uniform distribution of deaths

print(ax(udd, n=10, x=60, h=1/12, delta=log(1.05)))
print(Ax(udd, n=10, x=60, h=1/12, delta=log(1.05)))

tmp1 <- 1/(60*log(1.05)) * (1 - exp(-log(1.05)*(10)))
print(tmp1)

# Disability model

print(ax(dm, n=10, x=60, h = 1/12, delta=log(1.05)))
print(Ax(dm, n=10, x=60, h = 1/12, delta=log(1.05)))
