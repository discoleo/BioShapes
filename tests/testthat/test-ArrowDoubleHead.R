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

##### Double ArrowHead #####

cat("\nStarting: Double-Lined ArrowHead\n")

###### default join: join through #####
# - comparison between different dV (d.head argument);

###### Test 1 ######
x = c(0, 6); y = c(1, 6);
d = -1; dV = c(-1.5,1.5);
plot.base()
a1 = arrowDouble(x, y, d=d, dV=dV, lwd=2);
a2 = arrowDouble(c(x[1], 5), c(y[1], y[1]), d=d, dV=dV, lwd=2);
a3 = arrowDouble(c(x[1], x[1]), c(y[1], 5), d=d, dV=dV, lwd=2);
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
testArrow(h=h1, d=d, dV=dV)
testArrow(h=h2, d=d, dV=dV)
testArrow(h=h3, d=d, dV=dV)
testArrow(h=h4, d=d, dV=dV)
testArrow(h=h5, d=d, dV=dV)
testArrow(h=h6, d=d, dV=dV)


###### Test 2 ######
x = c(0, 6); y = c(1, 6) + 1;
d=-0.35; dV = c(-1.5, 1.5);
d.head = -1.5;
plot.base()
a1 = arrowDouble(x, y, d=d, d.head=d.head, dV=dV, lwd=2);
a2 = arrowDouble(c(x[1], 8), c(y[1], y[1]), d=d, dV=dV, d.head=d.head, lwd=2);
a3 = arrowDouble(c(x[1], x[1]), c(y[1], 8), d=d, dV=dV, d.head=d.head, lwd=2);
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
testArrow(h=h1, d=d.head, dV=dV)
testArrow(h=h2, d=d.head, dV=dV)
testArrow(h=h3, d=d.head, dV=dV)


###### Test 3 ######
x = c(0, 6); y = c(1, 80);
d = -3; d.head = -2;
scale = (100/12)*aspect_ratio_max
plot.base(ylim = c(0,100))
a1 = arrowDouble(x, y, d=d, d.head=d.head, lwd=2, scale=scale);
a2 = arrowDouble(c(x[1], 5), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, scale=scale);
a3 = arrowDouble(c(x[1], x[1]), c(y[1], 50), d=d, d.head=d.head, lwd=2, scale=scale);
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
d = -1; d2 = -1.5;
plot.base()
a1 = arrowDouble(x, y, d=d, lwd=2, join=1);
a2 = arrowDouble(c(x[1], 5), c(y[1], y[1]), d=d, lwd=2, join=1);
a3 = arrowDouble(c(x[1], x[1]), c(y[1], 5), d=d, lwd=2, join=1);
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
testArrow(h=h1, d=d)
testArrow(h=h2, d=d)
testArrow(h=h3, d=d)


###### Test 2 ######
x = c(0, 6); y = c(1, 6) + 1;
d=-1.5; d.head = -1.5; dV = c(-1.5,1.5);
plot.base()
a1 = arrowDouble(x, y, d=d, d.head=d.head, dV=dV, lwd=2, join=1);
a2 = arrowDouble(c(x[1], 8), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, join=1);
a3 = arrowDouble(c(x[1], x[1]), c(y[1], 8), d=d, d.head=d.head, lwd=2, join=1);
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
testArrow(h=h1, d=d)
testArrow(h=h2, d=d)
testArrow(h=h3, d=d)


###### Test 3 ######
x = c(0, 6); y = c(1, 80);
d = -2; d.head = -3;
scale = (100/12)*aspect_ratio_max
plot.base(ylim = c(0,100))
a1 = arrowDouble(x, y, d=d, d.head=d.head, lwd=2, scale=scale, join=1);
a2 = arrowDouble(c(x[1], 5), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, scale=scale, join=1);
a3 = arrowDouble(c(x[1], x[1]), c(y[1], 50), d=d, d.head=d.head, lwd=2, scale=scale, join=1);
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
