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

cat("Starting: X Shape ArrowHead\n")

##### X Shape ArrowHead #####
x = c(0, 6); y = c(1, 6);
plot.base()
arrowX(x, y, d=-1, lwd=2);
arrowX(c(x[1], 5), c(y[1], y[1]), d=-1, lwd=2);
arrowX(c(x[1], x[1]), c(y[1], 5), d=-1, lwd=2);

###### Test 1 ######
x = c(0, 6); y = c(1, 6);
d = -1;
plot.base()
a1 = arrowX(x, y, d=d, lwd=2);
a2 = arrowX(c(x[1], 5), c(y[1], y[1]), d=d, lwd=2);
a3 = arrowX(c(x[1], x[1]), c(y[1], 5), d=d, lwd=2);
# Head
h1 = a1$Head[[1]]

h2 = a2$Head[[1]]

h3 = a3$Head[[1]]
# - visual aids:
linesAid(h1, h2, h3)
linesAid(h1, h2, h3, id=c(2,4))
linesAid(h1, h2, h3, id=c(1,4))
linesAid(h1, h2, h3, id=c(2,3))

# Total length = (d^2 + dV[1]^2) + (d^2 + dV[2]^2)
cat("Test 1\n")
testArrow(h=h1, d=d)
testArrow(h=h2, d=d)
testArrow(h=h3, d=d)


###### Test 2 ######
x = c(0, 6); y = c(1, 6) + 1;
d = -1.5; d.head = c(-d-1, d+1);
plot.base()
a1 = arrowX(x, y, d=d, lwd=2, d.head=d.head);
a2 = arrowX(c(x[1], 8), c(y[1], y[1]), d=d, d.head=d.head, lwd=2);
a3 = arrowX(c(x[1], x[1]), c(y[1], 8), d=d, d.head=d.head, lwd=2);
# Head
h1 = a1$Head[[1]]

h2 = a2$Head[[1]]

h3 = a3$Head[[1]]
# - visual aids:
linesAid(h1, h2, h3)
linesAid(h1, h2, h3, id=c(2,4))
linesAid(h1, h2, h3, id=c(1,4))
linesAid(h1, h2, h3, id=c(2,3))

cat("Test 2\n")

testArrow(h=h1, d=d, dV=d.head)
testArrow(h=h2, d=d, dV=d.head)
testArrow(h=h3, d=d, dV=d.head)


###### Test 3 ######
x = c(0, 6); y = c(1, 80);
d = -2; d.head = c(-d-3, d+3);
scale = (100/12)*aspect_ratio_max
plot.base(ylim = c(-50,100))
a1 = arrowX(x, y, d=d, d.head=d.head, lwd=2, scale=scale);
a2 = arrowX(c(x[1], 5), c(y[1], y[1]), d=d, d.head=d.head, lwd=2, scale=scale);
a3 = arrowX(c(x[1], x[1]), c(y[1], 50), d=d, d.head=d.head, lwd=2, scale=scale);
# Head
h1 = a1$Head[[1]]

h2 = a2$Head[[1]]

h3 = a3$Head[[1]]
# - visual aids:
linesAid(h1, h2, h3)
linesAid(h1, h2, h3, id=c(2,4))
linesAid(h1, h2, h3, id=c(1,4))
linesAid(h1, h2, h3, id=c(2,3))

cat("Test 3: only visual\n")


# X: Different lengths
x = c(0, 6); y = c(1, 8);
plot.base()
arrowX(x, y, d=c(0.75, 2.0), d.head=c(-1,1), lwd=2)
abline(v=x[2], col="green");
points(x[2], y[2], col="red")
