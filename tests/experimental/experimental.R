###################
#
# Bachelor Thesis
#
# Title: BioShapes
#
# Candidate: Adrian Cotoc
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#
# in collaboration with Syonic SRL
# continues the work of Darian Voda
#
# GitHub: https://github.com/Adi131313/BioShapes

# TODO: create complex examples in R/Examples.R
# TODO: Centers for muscles
# TODO: Genome for examples.virus


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

################

### DNA Helix

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

h2 = dna.new(c(1,8), c(1,6), phi = c(0, pi), n.lin=6)
plot.base()
lines(h2, lwd=1)

h2 = dna.new(c(1,8), c(1,6), phi = c(pi/2, pi), n.lin=5)
plot.base()
lines(h2, lwd=2)

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
phi = c(pi, -pi/2); n = 3;
example.dna.tests(n=n, phi=phi)

###
phi = c(pi, -pi/2); n = 3.4;
example.dna.tests(n=n, phi=phi)

###
phi = c(pi, -pi/2); n = 3;
example.dna.tests(n=n, phi=phi, y0 = c(0, -2))

###
phi = c(pi, -pi/2); n = 1.7;
example.dna.tests(n=n, phi=phi, x0 = c(0, 6), y0 = c(0, -2))
# Note: (x0, y0) are NOT exact;
abline(v=6, col="blue")


### TODO: Design of monocyte
R = 20; n = 10;
lim = c(-R, R) + c(-2,2);
plot.base(xlim = lim, ylim = lim)
abline(h=0, col="green", lty=2)
xy = helix.rad(R = R, n=n, phi=pi, r = 2);
lines(xy, col="red")


#########################

### Neuron

# Test function
# TODO: Also one for Virus
description.neuron()

### Synapse
plot.base()
tmp = neuron(c(0,2), phi=pi/6, type.syn = "D")
lines(tmp)

###
plot.base()
tmp = neuron(c(0,2), phi=pi/6)
lines(tmp)

###
plot.base()
tmp = neuron(c(0,2), phi=pi/6, type.syn = "Tree")
lines(tmp)
# TODO: BUG for tree
lines(synapse(c(tmp[[6]]$x[2], tmp[[6]]$y[2]), slope=tan(pi/6), type="Tree", l=1.8, alpha=120))

### Test Neuron Synapses
plot.neuron = function(p, phi, type, dphi = c(0, -pi/9, pi/9)) {
  for(dp in dphi) {
    tmp = neuron(p, phi= phi + dp, type=type);
    lines(tmp);
  }
}

type = "Tree"
plot.base(ylim = c(-10, 10))
plot.neuron(c(4,8), phi = -pi/2, type=type)
plot.neuron(c(8,4), phi = pi, type=type)
plot.neuron(c(4,-8), phi = pi/2, type=type)
plot.neuron(c(0,-4), phi = 0, type=type)

### TODO 1: Dendrite as tree

######################

### Test: Bug in Function helix
plot.base();
p1 = c(3,1); p2 = c(1,5)
xy = helix(p1, p2); lines(xy);
points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red")
#
p1 = c(4,6); p2 = c(1, 5)
xy = helix(p1, p2); lines(xy);
points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red")
#
p1 = c(4,6); p2 = c(6,1)
xy = helix(p1, p2); lines(xy);
points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red")
#
p2 = c(4,6); p1 = c(1, 5)
xy = helix(p1, p2); lines(xy);
points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red")
