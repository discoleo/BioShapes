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


#######################

### Polycyclic Polyenes

lim = c(-5,5)
plot.base(xlim=lim, ylim=lim, axt=NULL, asp=1)
lines(polycycle.polyene(14, 7))
lines(annotate.polycycle(14))


################

################
### Analysis ###

### Newman Projection

###
test.proj.newman()

###
test.proj.newman(phi = c(0, pi/3), ligand = "CH3,H3C,F|Me,t-Bu,HO")


# compact
plot.base()
lines(proj.newman(c("F,H,H|OH,CH3,H"), pos = list(c(4,3,1), c(3,2,1))))

# explicit
plot.base()
lines(proj.newman(c("F", "H", "H"), c("OH", "CH3", "H"), pos = list(c(4,3,1), c(3,2,1))))

