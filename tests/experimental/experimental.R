#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. BSc Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#   GitHub: https://github.com/Adi131313/BioShapes
#
# 2. BSc Thesis: Darian Voda (2021-2022)



# TODO: create complex examples in R/Examples.R
# TODO: Centers for muscles
# TODO: Genome for example.virus
# TODO: description for Virus


######################

### Basic Functions

# Lines: Intersection
test.lines.simple(cbind(c(0,0,0,0), c(1,5,1,4)))
test.lines.simple(cbind(c(0,0,0,4), c(6,8,7,7)), add = TRUE)
test.lines.simple(cbind(c(1,1,1,7), c(1,5,1,4)), add = TRUE)
test.lines.simple(cbind(c(3,8,3,8), c(1,8,1,4)), add = TRUE)
test.lines.simple(cbind(c(2,2,2.5,4), c(2,6,4,6)), add = TRUE)
test.lines.simple(cbind(c(5,8,8,8), c(1,1,1,3)), add = TRUE)


######################

### Virus
N = 12
R = 3
lwd = 10
# Virus with ngons
virus = virus(N, R = R, center = c(5, 5), lwd = lwd)
plot.base()
lines(virus)
text(x = 5, y = 11, "Virus");
arrowDiamond(x=c(-2,1), y=c(4,4), lwd=2, d=-1)
arrowSimple(x=c(-1,2), y=c(-2,1), lwd=2, d=-1)

# Virus with circles
virus = virus(N, R = R, center = c(5, 5), ngon.spike = 0, lwd = lwd)
plot.base()
lines(virus)

### Virus 2 spikes

# Modified colors, lwd for virus and spikes
# length of spikes c(0, 2), R = 2
plot.base(xlim = c(-6, 6), ylim = c(-6, 6))
tmp = virus2()
lines(tmp)

### Virus & Genome

### DNA: ds Helix
cc = c(4,5); r = 3;
plot.base()
lines(virus(r, cc))
lst = genome(r * 0.75, cc)
lines(lst)

### Circle Arc
cc = c(4,5); r = 3;
plot.base()
lines(virus(r, cc))
lst = genome(r * 0.75, cc, type = "arc")
lines(lst)

### Circle
cc = c(4,5); r = 3;
plot.base()
lines(virus(r, cc))
lst = genome(r * 0.5, cc, type = "circle")
lines(lst)

### Simple Helix
cc = c(4,5); r = 3;
plot.base()
lines(virus2(r, cc))
lst = genome(r * 0.75, cc, type="helix", lwd=2)
lines(lst)

### css Helix
cc = c(4,5); r = 3;
plot.base()
lines(virus(r, cc))
lst = genome(r * 0.5, cc, type="css", lwd=2)
lines(lst)

### cds Helix / cyclic DNA
cc = c(4,5); r = 3;
plot.base()
lines(virus(r, cc))
lst = genome(r * 0.5, cc, type="cds", lwd=2)
lines(lst)


#################
### DNA Helix ###

# old algorithm
tmp1 = helix(c(0, 0), c(5, 5), n = 2, phi = 0, A = 3/4)
tmp2 = helix(c(0, 0), c(5, 5), n = 2, phi = pi, A = 3/4)
plot.base()
lines(tmp1, col = "red", lwd = 5)
lines(tmp2, col = "green", lwd = 5)

tmp = lapply(seq(40), function(id) {
  id = 3.5*id;
  lines(c(tmp1[[1]]$x[id], tmp2[[1]]$x[id]),
      c(tmp1[[1]]$y[id], tmp2[[1]]$y[id]), col = "blue");
})

pp = which.intersect.sin(c(0, pi), 2)
p0 = pp$x0 / (2*2*pi) * 5
abline(v = p0, col="purple")


### DNA

par.old = par(mfrow = c(2,1))

h2 = dna.new(c(1,8), c(1,6), phi = c(0, pi), n.lin=6)
plot.base(ylim = c(0, 8))
lines(h2)

h2 = dna.new(c(1,8), c(1,6), phi = c(pi/2, pi), n.lin=5, lwd=2)
plot.base(ylim = c(0, 8))
lines(h2)

par(par.old)


# TODO: Example showing design of DNA
# Intersection of 2 shifted-Sine Functions
# sin(x) = sin(x + phi)
phi = c(0.7, 1.5) * pi;
curve(sin(x + phi[1]), -2*pi, 6*pi, lwd=1.5)
curve(sin(x + phi[2]), add = T, col="red", lwd=1.5)
pp = which.intersect.sin(phi, 3)
abline(v = pp$x0, col="blue")
points(pp$x0, sin(pp$x1), col="green", lwd=2)


### Tests for DNA

###
phi = c(pi, -pi/2);
test.dna.nb(phi);

phi = c(pi/2, pi);
test.dna.nb(phi);


############

### TODO: Design of monocyte
R = 20; n = 10;
lim = c(-R, R) + c(-2,2);
plot.base(xlim = lim, ylim = lim)
abline(h=0, col="green", lty=2)
xy = helix.rad(R = R, n=n, phi=pi, r = 2);
lines(xy, col="red")


#########################

### Neuron

### Test Neuron Synapses
plot.neuron = function(p, phi, type, dphi = c(0, -pi/9, pi/9)) {
  for(dp in dphi) {
    tmp = neuron(p, phi = phi + dp, type=type);
    lines(tmp);
  }
}

###
par.old = par(mfrow = c(2,2))

### Synapse
plot.base()
tmp = neuron(c(0,2), phi=pi/6, type.syn = "D")
lines(tmp)

###
plot.base()
tmp = neuron(c(0,2), phi=pi/6)
lines(tmp)

###
tmp = neuron(c(0,2), phi=pi/6, type.syn = "Tree", synapse = list(l = 2))
plot.base()
lines(tmp)

###
type = "Tree"
plot.base(ylim = c(-10, 10))
plot.neuron(c(4,8), phi = -pi/2, type=type)
plot.neuron(c(8,4), phi = pi, type=type)
plot.neuron(c(4,-8), phi = pi/2, type=type)
plot.neuron(c(0,-4), phi = 0, type=type)

par(par.old);


### TODO: Dendrite as tree

######################

### Test: Bug in Function helix
test.helix.directions()


###################

### Tests: Pins

test.pins()


###################

### DH Arrows

test.arrow.Half()


### Tests: Boxes

test.box.cap()


### Box with Elliptic Cap

###
par.old = par(mfrow = c(2,2))

test.box.capEllipse()

test.box.capEllipse(fill = "red")

test.box.capEllipse(scale = 2)

test.box.capEllipse(scale = 1/2)

par(par.old)


###################

### Ig-like Domains

# minimal test:
plot.base()
tt = c(1/2, 3/4); l = 1/2;
lines(mol.IgDomain(c(1,1), c(1,8), t = tt, l=l))
lines(mol.IgDomain(c(2,3), c(1,8), t = tt, l=l, lwd=2, col.ig="red"))
lines(mol.IgDomain(c(3,5), c(1,8), t = tt, l=l, lwd=3, col.ig=c("blue", "red")))


### Ig Monomer

# Upper quadrants:
plot.base(xlim = c(-10, 10))
lines(mol.Ig(c( 0,1)))
lines(mol.Ig(c(-2,1), phi = pi - pi/10))
lines(mol.Ig(c( 2,1), phi = pi/10))

# Lower Quadrants:
plot.base(xlim = c(-10, 10))
lines(mol.Ig(c(-3,6), phi = pi + pi/12))
lines(mol.Ig(c(-1,3), phi = pi + pi/4))
lines(mol.Ig(c( 0,1), phi = 3*pi/2))
lines(mol.Ig(c( 1,3), phi = 2*pi - pi/4))
lines(mol.Ig(c( 3,6), phi = 2*pi - pi/12))


######################

### Lines / Connectors

par.old = par(mfrow = c(2,2))

plot.base()
x = c(1,6); y = c(1,6)
lines(connect.lines.stairs(x, y, slope=1/2))
points(x, y, col="red")


plot.base()
x = c(0, 8); y = c(0, 8)
tmp = sapply(c(1.2, 1.6, 2.2, 3.5), \(id) {
	lines(connect.lines.stairs(x, y, slope = id, t = 1/id))
})
points(x, y, col="red")


###
plot.base()
x = c(1,6); y = c(1,6)
lines(connect.lines.stairsEq(x, y))
points(x, y, col="red")

x = c(-1, 4, 9, 4, -1); y = c(4, 9, 4, -1, 4);
plot.base()
for(i in seq(4)) {
	id = c(i, i+1);
	lines(connect.lines.stairsEq(x[id], y[id]))
	lines(connect.lines.stairsEq(x[id], y[id], up = FALSE))
}
points(x, y, col="red")


par(par.old)

