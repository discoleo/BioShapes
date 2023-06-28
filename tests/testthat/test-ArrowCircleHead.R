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

cat("Starting: Circle ArrowHead\n")

##### Circle ArrowHead #####
x = c(0, 6); y = c(1, 6);
plot.base()
arrowCircle(x, y, r=0.5, lwd=2);
arrowCircle(c(x[1], 5), c(y[1], y[1]), r=0.5, lwd=2);
arrowCircle(c(x[1], x[1]), c(y[1], 5), r=0.5, lwd=2);

###### Test 1 ######
x = c(0, 6); y = c(1, 6);
r=0.5
plot.base()
a1 = arrowCircle(x, y, r=r, lwd=2);
a2 = arrowCircle(c(x[1], 5), c(y[1], y[1]), r=r, lwd=2);
a3 = arrowCircle(c(x[1], x[1]), c(y[1], 5), r=r, lwd=2);
# Head
h1 = a1$Head[[1]]

h2 = a2$Head[[1]]

h3 = a3$Head[[1]]
# - visual aids:
points(x[2], y[2], col="green", lwd=2)
points(x[1], 5, col="green", lwd=2)
points(5, y[1], col="green", lwd=2)
# Total length = (d^2 + dV[1]^2) + (d^2 + dV[2]^2)
cat("Test 1\n")
testArrow(h=h1, d=d)
testArrow(h=h2, d=d)
testArrow(h=h3, d=d)


###### Test 2 ######
x = c(0, 6); y = c(1, 6) + 1;
r = 1.5;
plot.base()
a1 = arrowCircle(x, y, r=r, lwd=2, join=3);
a2 = arrowCircle(c(x[1], 8), c(y[1], y[1]), r=r, lwd=2, join=2);
a3 = arrowCircle(c(x[1], x[1]), c(y[1], 8), r=r, lwd=2);
# Head
h1 = a1$Head[[1]]

h2 = a2$Head[[1]]

h3 = a3$Head[[1]]
# - visual aids:
linesAid(h1, h2, h3)

cat("Test 2\n")

testArrow(h=h1, d=d, dV=d.head)
testArrow(h=h2, d=d, dV=d.head)
testArrow(h=h3, d=d, dV=d.head)


###### Test 3 ######
x = c(0, 6); y = c(1, 80);
r = 2;
scale = (100/12)*aspect_ratio_max
plot.base(ylim = c(-50,100))
a1 = arrowCircle(x, y, r=r, lwd=2, scale=scale);
a2 = arrowCircle(c(x[1], 5), c(y[1], y[1]), r=r, lwd=2, scale=scale);
a3 = arrowCircle(c(x[1], x[1]), c(y[1], 50), r=r,lwd=2, scale=scale);
# Head
h1 = a1$Head[[1]]

h2 = a2$Head[[1]]

h3 = a3$Head[[1]]
# - visual aids:
linesAid(h1, h2, h3)

cat("Test 3: only visual\n")

###### Filled Circles ######
x = c(0, 6); y = c(1, 6);
r = 0.5;
fill = "#E0B0B0";
plot.base()
arrowCircle(x, y, r=r, lwd=2, fill=fill);
arrowCircle(c(x[1], 5), c(y[1], y[1]), r=r, lwd=2, fill=fill);
arrowCircle(c(x[1], x[1]), c(y[1], 5), r=r, lwd=2, fill=fill);
