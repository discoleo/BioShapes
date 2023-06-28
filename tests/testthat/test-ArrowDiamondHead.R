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

cat("Starting: Diamond ArrowHead\n")

##### Diamond ArrowHead #####

# default join
x = c(0, 6); y = c(1, 6);
plot.base()
arrowDiamond(x, y, d=-1, lwd=2);
arrowDiamond(c(x[1], 5), c(y[1], y[1]), d=-1, lwd=2);
arrowDiamond(c(x[1], x[1]), c(y[1], 5), d=-1, lwd=2);

# join through
x = c(0, 6); y = c(1, 6);
plot.base()
arrowDiamond(x, y, d=-1, lwd=2, join=2);
arrowDiamond(c(x[1], 5), c(y[1], y[1]), d=-1, lwd=2, join=2);
arrowDiamond(c(x[1], x[1]), c(y[1], 5), d=-1, lwd=2, join=2);


###### Test 1 ######
x = c(0, 6); y = c(1, 6);
d = -1;
plot.base()
a1 = arrowDiamond(x, y, d=d, lwd=2);
a2 = arrowDiamond(c(x[1], 5), c(y[1], y[1]), d=d, lwd=2);
a3 = arrowDiamond(c(x[1], x[1]), c(y[1], 5), d=d, lwd=2);
# Head
h1 = a1$Head[[1]]
h2 = a2$Head[[1]]
h3 = a3$Head[[1]]
# - visual aids:
linesAid(h1, h2, h3)

cat("Test 1\n")
testArrow(h=h1, d=d)
testArrow(h=h2, d=d)
testArrow(h=h3, d=d)


###### Test 2 ######
# - narrow ArrowHead;
x = c(0, 6); y = c(1, 6) + 1;
d=-0.5; d.head = c(-1.5, 1.5);
plot.base()
a1 = arrowDiamond(x, y, d=d, d.head=d.head, lwd=2);
a2 = arrowDiamond(c(x[1], 5), c(y[1], y[1]), d=d, d.head=d.head, lwd=2);
a3 = arrowDiamond(c(x[1], x[1]), c(y[1], 8), d=d, d.head=d.head, lwd=2);
# Head
h1 = a1$Head[[1]]
h2 = a2$Head[[1]]
h3 = a3$Head[[1]]
# - visual aids:
linesAid(h1, h2, h3)

cat("Test 2\n")
# Total length = (d^2 + dV[1]^2) + (d^2 + dV[2]^2)
testArrow(h=h1, d=d, dV=d.head)
testArrow(h=h2, d=d, dV=d.head)
testArrow(h=h3, d=d, dV=d.head)


###### Test 3 ######
x = c(0, 4); y = c(1, 60);
x2 = c(0, 6); y2 = c(1, 20);
d = -2; d.head = c(-1/2, 1/2);
scale = (100/12) * aspect_ratio_max;
plot.base(ylim = c(0,100))
a1 = arrowDiamond(x, y, d=d, d.head=d.head, lwd=2, scale=scale);
a2 = arrowDiamond(x2, y2, d=d, d.head=d.head, lwd=2, scale=scale);
a3 = arrowDiamond(c(x[1], 9), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, scale=scale);
a4 = arrowDiamond(c(x[1], x[1]), c(y[1], 50), d=d, d.head=d.head, lwd=2, scale=scale);
# Head
h1 = a1$Head[[1]]
h2 = a2$Head[[1]]
h3 = a3$Head[[1]]
h4 = a4$Head[[1]]
# - visual aids:
linesAid(h1, h2, h3, h4)
cat("Test 3: only visual\n\n")

