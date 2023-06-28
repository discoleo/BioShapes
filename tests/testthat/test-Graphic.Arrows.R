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

#### Multiple-lined ArrowHead ####
# TODO: Junction
plot.base()
x = c(0, 6); y = c(1, 6);
arrowN(x, y, n = 5, d=0.5, lwd=2)

#### Double Lined Inverted ArrowHead ####
plot.base()
x = c(0, 6); y = c(1, 6);
arrowDoubleInverted(x, y, lwd=2);
abline(v=x[2], col="green");
points(x[2], y[2], col="red")


#### Double lined ArrowHead ####
# TODO: Junction
plot.base()
x = c(0, 6); y = c(1, 6);
arrowDouble(x, y, d=0.5, lwd=2);
abline(v=x[2], col="green");
points(x[2], y[2], col="black")

#### Circle Arrowhead ####
plot.base()
x = c(0, 6); y = c(1, 6);
arrowCircle(x, y, r=0.5, lwd=2);
abline(v=x[2], col="green");
points(x[2], y[2], col="red")

#### Solid Square ArrowHead ####
plot.base()
x = c(0, 6); y = c(1, 6);
arrowSolidSquare(x, y, d=-2, lwd=2);
abline(v=x[2], col="green");
points(x[2], y[2], col="red")

#### Triangle  ArrowHead ###

### Test 1
x = c(0, 6); y = c(1, 6);
d = -1;
plot.base()
a1 = arrowTriangle(x, y, d=d, lwd=2);
a2 = arrowTriangle(c(x[1], 5), c(y[1], y[1]), d=d, lwd=2);
a3 = arrowTriangle(c(x[1], x[1]), c(y[1], 5), d=d, lwd=2);
# Head
h1 = a1$Head[[1]]
h2 = a2$Head[[1]]
h3 = a3$Head[[1]]
# 4 edges: d * c(sqrt(2), sqrt(2), 0, 2);
# => edge 1 is counted twice; edge 4 is not counted;
stopifnot(round(Dsquare(h1, h1$x[2], h1$y[2]) - 6*d^2, 8) == 0)
stopifnot(round(Dsquare(h2, h2$x[2], h2$y[2]) - 6*d^2, 8) == 0)
stopifnot(round(Dsquare(h3, h3$x[2], h3$y[2]) - 6*d^2, 8) == 0)


###########################

#### Square-Wave Arrow ####
x = c(1,5); y = c(1,7);
plot.base()
arrowSquareWave(x, y, n=5)
arrowSquareWave(x + c(0, 4), y - c(0, 3), n=6, col="red")

# No. of Teeth: Variations
x = c(-1, 8); y = c(0, 0);
dy = 2.5;
plot.base()
arrowSquareWave(x, y, n=5)
arrowSquareWave(x, y + 1*dy, n=6)
arrowSquareWave(x, y + 2*dy, n=7)
arrowSquareWave(x, y + 3*dy, n=8)


