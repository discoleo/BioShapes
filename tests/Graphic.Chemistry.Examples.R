#### Poly-cyclic chemical molecules ####


par.old = par(mfrow = c(3,2))

### Cholesterol backbone
plot.base(ylim = c(-2,4))
lines(parseCycles("6|6\\6|5"))

### [1] Pentacyclic Triterpene backbone
plot.base(ylim = c(-2,4))
lines(parseCycles("6|6\\6|6\\6"))

### [3] Molecules from vegetables
### Psoralens & Pabulenol quasi-backbone
# - vegetables: parsnip (pastarnac) & parsley (patrunjel);
plot.base(ylim = c(-2,4))
lines(parseCycles("6|6|5"))

### Spiro-Derivatives
plot.base(ylim = c(-2,4))
lines(parseCycles("6|6|5<5"))

plot.base(ylim = c(-2,4))
lines(parseCycles("6|6|5<7"))

### Other
plot.base(ylim = c(-2,4))
lines(parseCycles("6|6|6|5<3"))

par(par.old)

###
plot.base()
lines(parseCycles("5-5"))

### Sesquifulvalene backbone
plot.base()
lines(parseCycles("7=5"))


################

################
### Analysis ###

### Newman Projection

plot.base()
lines(proj.newman(c("F", "H", "H"), c("OH", "CH3", "H"), pos = list(c(4,3,1), c(3,2,1))))

