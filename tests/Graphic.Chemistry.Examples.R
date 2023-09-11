#### Poly-cyclic chemical molecules ####

#### Cholesterol backbone ####
plot.base()
lines(parseCycles("6|6\\6|5"))

##### [1] Pentacyclic Triterpene backbone ####
plot.base()
lines(parseCycles("6|6\\6|6\\6"))

#### [3] Molecules from vegetables; ####
##### Psoralens & Pabulenol quasi-backbone #####
# - vegetables: parsnip (pastarnac) & parsley (patrunjel);
plot.base()
lines(parseCycles("6|6|5"))

#### Spiro-Derivatives ####
plot.base()
lines.arrow(list(parseCycles("6|6|5<5"), list()))

plot.base()
lines.arrow(list(parseCycles("6|6|5<7"), list()))

#### Other ####
plot.base()
lines.arrow(list(parseCycles("6|6|6|5<3"), list()))

plot.base()
lines.arrow(list(parseCycles("5-5"), list()))

#### Sesquifulvalene backbone ####
plot.base()
lines.arrow(list(parseCycles("7=5"), list()))

################

################
### Analysis ###

### Newman Projection

plot.base()
lines(proj.newman(c("F", "H", "H"), c("OH", "CH3", "H"), pos = list(c(4,3,1), c(3,2,1))))

