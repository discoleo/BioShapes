###################
#
# Thesis
#
# Title: BioShapes
#
# Candidate: Darian Voda
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#
# in collaboration with Syonic SRL
#
# GitHub: https://github.com/DarianFlorianVoda/Diagram-Generator


###############
#### Tests ####

#### Diagonal, H & V Tests ####

cat("Starting: Double Inverted ArrowHead\n")

##### Double Inverted ArrowHead #####
x = c(0, 6); y = c(1, 6);
plot.base()
arrowDoubleInverted(x, y, d=-1, lwd=2);
arrowDoubleInverted(c(x[1], 5), c(y[1], y[1]), d=-1, lwd=2);
arrowDoubleInverted(c(x[1], x[1]), c(y[1], 5), d=-1, lwd=2);

###### Test 1 ######
x = c(0, 6); y = c(1, 6);
d = -1; dH = 0.5; # default dH
plot.base()
a1 = arrowDoubleInverted(x, y, d=d, lwd=2);
a2 = arrowDoubleInverted(c(x[1], 5), c(y[1], y[1]), d=d, lwd=2);
a3 = arrowDoubleInverted(c(x[1], x[1]), c(y[1], 5), d=d, lwd=2);
# Head
h1 = a1$Head[[1]]
h2 = a1$Head[[2]]

h3 = a2$Head[[1]]
h4 = a2$Head[[2]]

h5 = a3$Head[[1]]
h6 = a3$Head[[2]]
# - visual aids:
linesAid(h1, h2, h3, h4, h5, h6)

cat("Test 1\n")
# Total length = (d^2 + dV[1]^2) + (d^2 + dV[2]^2)
testArrow(h=h1, d=dH)
testArrow(h=h2, d=dH)
testArrow(h=h3, d=dH)
testArrow(h=h4, d=dH)
testArrow(h=h5, d=dH)
testArrow(h=h6, d=dH)


###### Test 2 ######
x = c(0, 6); y = c(1, 6) + 1;
d=-1.75; dH=2.5
d.head = c(-dH+0.5,dH-0.5);
plot.base()
a1 = arrowDoubleInverted(x, y, d=d, dH=dH, d.head=d.head, lwd=2);
a2 = arrowDoubleInverted(c(x[1], 8), c(y[1], y[1]), d=d, dH=dH, d.head=d.head, lwd=2);
a3 = arrowDoubleInverted(c(x[1], x[1]), c(y[1], 8), d=d, dH=dH, d.head=d.head, lwd=2);
# Head
h1 = a1$Head[[1]]
h2 = a1$Head[[2]]

h3 = a2$Head[[1]]
h4 = a2$Head[[2]]

h5 = a3$Head[[1]]
h6 = a3$Head[[2]]
# - visual aids:
linesAid(h1, h2, h3, h4, h5, h6)

cat("Test 2\n")
testArrow(h=h1, d=dH, dV=d.head)
testArrow(h=h2, d=dH, dV=d.head)
testArrow(h=h3, d=dH, dV=d.head)


###### Test 3 ######
x = c(0, 6); y = c(1, 80);
d = -0.33; d.head = c(-2,2); dH = 0.5; # default dH
scale = (100/12)*aspect_ratio_max
plot.base(ylim = c(0,100))
a1 = arrowDoubleInverted(x, y, d=d, d.head=d.head, lwd=2, scale=scale);
a2 = arrowDoubleInverted(c(x[1], 5), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, scale=scale);
a3 = arrowDoubleInverted(c(x[1], x[1]), c(y[1], 50), d=d, d.head=d.head, lwd=2, scale=scale);
# Head
h1 = a1$Head[[1]]
h2 = a1$Head[[2]]

h3 = a2$Head[[1]]
h4 = a2$Head[[2]]

h5 = a3$Head[[1]]
h6 = a3$Head[[2]]
# - visual aids:
linesAid(h1, h2, h3, h4, h5, h6)

cat("Test 3: only visual\n\n")


###### join with innermost ">" #####
cat("Test: join\n")
###### Test 1 ######
x = c(0, 6); y = c(1, 6);
d = -1; dH = 0.5; # default dH
plot.base()
a1 = arrowDoubleInverted(x, y, d=d, lwd=2, join=2);
a2 = arrowDoubleInverted(c(x[1], 5), c(y[1], y[1]), d=d, lwd=2, join=2);
a3 = arrowDoubleInverted(c(x[1], x[1]), c(y[1], 5), d=d, lwd=2, join=2);
# Head
h1 = a1$Head[[1]]
h2 = a1$Head[[2]]

h3 = a2$Head[[1]]
h4 = a2$Head[[2]]

h5 = a3$Head[[1]]
h6 = a3$Head[[2]]
# - visual aids:
linesAid(h1, h2, h3, h4, h5, h6)
# Total length = (d^2 + dV[1]^2) + (d^2 + dV[2]^2)
cat("Test 1\n")
testArrow(h=h1, d=dH)
testArrow(h=h2, d=dH)
testArrow(h=h3, d=dH)


###### Test 2 ######
x = c(0, 6); y = c(1, 6) + 1;
d=-1.5; d.head = -1.5; dH = 0.5; # default dH
plot.base()
a1 = arrowDoubleInverted(x, y, d=d, d.head=d.head, lwd=2, join=2);
a2 = arrowDoubleInverted(c(x[1], 8), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, join=2);
a3 = arrowDoubleInverted(c(x[1], x[1]), c(y[1], 8), d=d, d.head=d.head, lwd=2, join=2);
# Head
h1 = a1$Head[[1]]
h2 = a1$Head[[2]]

h3 = a2$Head[[1]]
h4 = a2$Head[[2]]

h5 = a3$Head[[1]]
h6 = a3$Head[[2]]
# - visual aids:
linesAid(h1, h2, h3, h4, h5, h6)
cat("Test 2\n")
testArrow(h=h1, d=dH, dV=d.head)
testArrow(h=h2, d=dH, dV=d.head)
testArrow(h=h3, d=dH, dV=d.head)


###### Test 3 ######
x = c(0, 6); y = c(1, 80);
d = -2; d.head = -3;
scale = (100/12)*aspect_ratio_max
plot.base(ylim = c(0,100))
a1 = arrowDoubleInverted(x, y, d=d, d.head=d.head, lwd=2, scale=scale, join=2);
a2 = arrowDoubleInverted(c(x[1], 5), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, scale=scale, join=2);
a3 = arrowDoubleInverted(c(x[1], x[1]), c(y[1], 50), d=d, d.head=d.head, lwd=2, scale=scale, join=2);
# Head
h1 = a1$Head[[1]]
h2 = a1$Head[[2]]

h3 = a2$Head[[1]]
h4 = a2$Head[[2]]

h5 = a3$Head[[1]]
h6 = a3$Head[[2]]
# - visual aids:
linesAid(h1, h2, h3, h4, h5, h6)

cat("Test 3: only visual\n")
