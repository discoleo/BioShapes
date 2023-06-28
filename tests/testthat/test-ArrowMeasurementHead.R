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

cat("Starting: Measurement ArrowHead\n")

##### Measurement ArrowHead #####
x = c(0, 6); y = c(1, 6);
plot.base()
arrowMeasure(x, y, d=-1, lwd=2);
arrowMeasure(c(x[1], 5), c(y[1], y[1]), d=-1, lwd=2);
arrowMeasure(c(x[1], x[1]), c(y[1], 5), d=-1, lwd=2);

###### Test 1 ######
x = c(0, 6); y = c(1, 6);
d = -1;
plot.base()
a1 = arrowMeasure(x, y, d=d, lwd=2);
a2 = arrowMeasure(c(x[1], 5), c(y[1], y[1]), d=d, lwd=2);
a3 = arrowMeasure(c(x[1], x[1]), c(y[1], 5), d=d, lwd=2);
# Head
h1 = a1$Head[[1]]
h4 = a1$Head[[2]]
h5 = list()
h5$x = c(h4$x, h1$x)
h5$y = c(h4$y, h1$y)

h7 = list()
h2 = a2$Head[[1]]
h6 = a2$Head[[2]]
h7$x = c(h6$x, h2$x)
h7$y = c(h6$y, h2$y)

h9 = list()
h3 = a3$Head[[1]]
h8 = a3$Head[[2]]
h9$x = c(h8$x, h3$x)
h9$y = c(h8$y, h3$y)
# - visual aids:
linesAid(h1, h2, h3)
linesAid(h5, h7, h9, id=c(1,3))
linesAid(h5, h7, h9, id=c(2,5))
# Total length = (d^2 + dV[1]^2) + (d^2 + dV[2]^2)
cat("Test 1\n")
testArrow(h=h1, d=d)
testArrow(h=h2, d=d)
testArrow(h=h3, d=d)


###### Test 2 ######
x = c(0, 6); y = c(1, 6) + 1;
d = -1.5; #d.head = -3;
d.head = c(-d-2, d+2)
plot.base()
a1 = arrowMeasure(x, y, d=d, lwd=2, d.head=d.head);
a2 = arrowMeasure(c(x[1], 8), c(y[1], y[1]), d=d, d.head=d.head, lwd=2);
a3 = arrowMeasure(c(x[1], x[1]), c(y[1], 8), d=d, d.head=d.head, lwd=2);
# Head
h1 = a1$Head[[1]]
h4 = a1$Head[[2]]
h5 = c()
h5$x = c(h4$x, h1$x)
h5$y = c(h4$y, h1$y)

h7 = c()
h2 = a2$Head[[1]]
h6 = a2$Head[[2]]
h7$x = c(h6$x, h2$x)
h7$y = c(h6$y, h2$y)

h9 = c()
h3 = a3$Head[[1]]
h8 = a3$Head[[2]]
h9$x = c(h8$x, h3$x)
h9$y = c(h8$y, h3$y)
# - visual aids:
linesAid(h1, h2, h3)
linesAid(h5, h7, h9, id=c(1,3))
linesAid(h5, h7, h9, id=c(2,5))
cat("Test 2\n")
testArrow(h=h1, d=d, dV=d.head)
testArrow(h=h2, d=d, dV=d.head)
testArrow(h=h3, d=d, dV=d.head)


###### Test 3 ######
x = c(0, 6); y = c(1, 80);
d = -2; d.head = c(-d-3, d+3);
scale = (100/12)*aspect_ratio_max;
plot.base(ylim = c(-50,100))
a1 = arrowMeasure(x, y, d=d, d.head=d.head, lwd=2, scale=scale);
a2 = arrowMeasure(c(x[1], 5), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, scale=scale);
a3 = arrowMeasure(c(x[1], x[1]), c(y[1], 50), d=d, d.head=d.head, lwd=2, scale=scale);
# Head
h1 = a1$Head[[1]]
h4 = a1$Head[[2]]
h5 = c()
h5$x = c(h4$x, h1$x)
h5$y = c(h4$y, h1$y)

h7 = c()
h2 = a2$Head[[1]]
h6 = a2$Head[[2]]
h7$x = c(h6$x, h2$x)
h7$y = c(h6$y, h2$y)

h9 = c()
h3 = a3$Head[[1]]
h8 = a3$Head[[2]]
h9$x = c(h8$x, h3$x)
h9$y = c(h8$y, h3$y)
# - visual aids:
linesAid(h1, h2, h3)
linesAid(h5, h7, h9, id=c(1,3))
linesAid(h5, h7, h9, id=c(2,5))

cat("Test 3: only visual\n")

###### Measurement Arrow 1 line ######
x = c(0, 6); y = c(1, 6);
dT = c(-1); dV = c(-1,1) / 2;
plot.base()
arrowMeasure(x, y, d=-1, d.head=dV, dT=dT, lwd=2);
arrowMeasure(c(x[1], 5), c(y[1], y[1]), d=-1, d.head=dV, dT=dT, lwd=2);
arrowMeasure(c(x[1], x[1]), c(y[1], 5), d=-1, d.head=dV, dT=dT, lwd=2);
