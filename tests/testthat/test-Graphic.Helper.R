#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis: Darian Voda (2021-2022)
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#   [old] GitHub: https://github.com/DarianFlorianVoda/Diagram-Generator


# TODO:
# - test cases for various types of base-lines and points;
# - base line:
#  -- origin: (0, 0), (0, 1), (0, 2), (0, -2), ...;
#  -- slope: ascending X, descending X, vertical X, horizontal X,
#     slightly perturbed vertical or horizontal X;
# - separate source file with the various tests X;

# Base-Line
# TODO: various origins, e.g.(x0=0, y0=0) vs non-zero;
# Note: x = (x0, x1); y = (y0, y1);
# Testing various origins

### Test:

testReflection = function(p0, x, y) {
  p = reflect(x, y, p0)
  points(p0[1], p0[2])
  points(p[1], p[2]);
  lines(c(p0[1], p[1]), c(p0[2], p[2]), col="green")
  # TODO: stopifnot(p == p0);
  return(p); # others: invisible(p);
}

testShiftPoint = function(p0, x, y, d=1, color="green") {
  p = shiftPoint(p0, x, y, d)
  print(p)
  points(p0[1], p0[2])
  points(p, col=color)
}


testShiftLine = function(x, y, d=1, color = "orange") {
  l = shiftLine(x, y, d)
  len = nrow(l) / 2;
  sapply(seq(len), function(id) lines(l[c(id, id+len),1], l[c(id, id+len),2], col=color))
}

testRoundValues = function(x, y, p, decimal=6){
  div = (x[2]-x[1])**2 + (y[2]-y[1])**2
  value = p*div
  paste("c(", paste(value, "/", div, sep="", collapse=", "), "),", decimal)
}

##### Test: Reflections #####

###### Horizontal Reflection ######

# Line Start = (0, 4), End = (5, 4)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(4, 4);
lines(x, y, lwd=2, col="red")
abline(v=x, col="green", lty=3)

# points(x, y)!
# Test 1
p0 = c(2, 6)
p = testReflection(p0, x, y)
stopifnot(p==c(2, 2))

# Test 2
p0 = c(2.5, 2)
p = testReflection(p0, x, y)
stopifnot(p==c(2.5, 6))

# Test 3
p0 = c(1, 1)
p = testReflection(p0, x, y)
stopifnot(p == c(1, 7))

# Test 4
p0 = c(4, 9)
p = testReflection(p0, x, y)
stopifnot(p == c(4, -1))

# Test 5
p0 = c(6, 6)
p = testReflection(p0, x, y)
stopifnot(p==c(6, 2))

# Test 6
p0 = c(-1, -1)
p = testReflection(p0, x, y)
stopifnot(p==c(-1, 9))

# Test 7
p0 = c(3, 4)
p = testReflection(p0, x, y)
stopifnot(p==c(3, 4))

# Test 8
p0 = c(3.5, 4.25)
p = testReflection(p0, x, y)
stopifnot(p==c(3.5, 3.75))

cat("Finished Subsection [Reflections]: Horizontal Reflection - part 1\n");

###### Origin: (0, y) ######

# Tests: (0, 0), (0, 1), (0, 2), (0, -2)

# Line Start = (0, 0), End = (5, 4)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(0, 4);
lines(x, y, lwd=2, col="red")
abline(v=x, col="green", lty=3)

# points(x, y)!
# Test 1
p0 = c(0, 0)
p = testReflection(p0, x, y)
stopifnot(p == c(0, 0))

# Test 2
p0 = c(4, 0)
p = testReflection(p0, x, y)
stopifnot(round(p - c(36 / 41, 3 + 37/41), 12) == 0)

# Test 3
p0 = c(0, 2)
p = testReflection(p0, x, y)
stopifnot(round(p - c(80/41, - 18/41), 12) == 0)

# Test 4
p0 = c(3.75, 0.3)
p = testReflection(p0, x, y)
stopifnot(round(p - c( 45.75/41, 147.3/41 ), 6) == 0)

# Line Start = (0, 1), End = (0, -2)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0,0); y = c(1, -2);
lines(x, y, lwd=2, col="red")

# Test 1
p0 = c(-1, -1)
p = testReflection(p0, x, y)
stopifnot(p == c(1, -1))

# Test 2
p0 = c(0, 0)
p = testReflection(p0, x, y)
stopifnot(p == c(0, 0))

cat("Finished Subsection [Reflections]: Origin: (0, y) - part 2\n");

###### Vertical Reflection ######

# Line Start = (3, 1), End = (3, 8)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(3, 3); y = c(1, 8);
lines(x, y, lwd=2, col="red")
abline(h=y, col="green", lty=3)

# points(x, y)!
# Test 1
p0 = c(4, 7)
p = testReflection(p0, x, y)
stopifnot(all(p == c(2,7)))

# Test 2
p0 = c(6, 5)
p = testReflection(p0, x, y)
stopifnot(all(p == c(0,5)))

# Test 3
p0 = c(8, 3)
p = testReflection(p0, x, y)
stopifnot(all(p == c(-2,3)))

# Test 4
p0 = c(2, 6)
p = testReflection(p0, x, y)
stopifnot(all(p == c(4,6)))

# Test 5
p0 = c(3, 4)
p = testReflection(p0, x, y)
stopifnot(all(p == c(3,4)))

# Test 6
p0 = c(2, 0)
p = testReflection(p0, x, y)
stopifnot(all(p == c(4,0)))

cat("Finished Subsection [Reflections]: Vertical Reflection - part 3\n");

###### Ascending Reflection #######

# Line Start = (0, 0), End = (6, 6)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0,6); y = c(0, 6);
lines(x, y, lwd=2, col="red")

# points(x, y)!
# Test 1
p0 = c(1, 2)
p = testReflection(p0, x, y)
stopifnot(all(p == c(2,1)))

# Test 2
p0 = c(2, 4)
p = testReflection(p0, x, y)
stopifnot(all(p == c(4,2)))

# Test 3
p0 = c(4, 7)
p = testReflection(p0, x, y)
stopifnot(all(p == c(7,4)))

# Test 4
p0 = c(1, 0)
p = testReflection(p0, x, y)
stopifnot(all(p == c(0,1)))

# Test 5
p0 = c(4, 4)
p = testReflection(p0, x, y)
stopifnot(all(p == c(4,4)))


# Test 6
p0 = c(0, -1)
p = testReflection(p0, x, y)
stopifnot(all(p == c(-1,0)))

cat("Finished Subsection [Reflections]: Ascending Reflection - part 4\n");

###### Descending Reflection ######

# Line Start = (0, 8), End = (8, 1)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 8); y = c(8, 1);
lines(x, y, lwd=2, col="red")
abline(v=x, col="green", lty=3)

# points(x, y)!
# Test 1
p0 = c(1, 8)
p = testReflection(p0, x, y)
stopifnot(round(p - c(0, 7) - c(15,1)/113, 12) == 0)

# Test 2
p0 = c(3, 8)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 45/113, 568/113 ), 6) == 0)

# Test 3
p0 = c(6, 7)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 202/113, 247/113 ), 6) == 0)

# Test 4
p0 = c(4, 2)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 732/113, 546/113 ), 6) == 0)

# Test 5
p0 = c(8, 1)
p = testReflection(p0, x, y)
stopifnot((all(p == c(8,1))))

cat("Finished Subsection [Reflections]: Descending Reflection - part 5\n");

###### Slightly Perturbed Vertical Reflection ######

# Line Start = (3.25, 1), End = (3, 8)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(3.25, 3); y = c(1, 8);
lines(x, y, lwd=2, col="red")
abline(h=y, col="green", lty=3)

# points(x, y)!
# Test 1
p0 = c(4, 7)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 101.75/49.0625, 340.0625/49.0625 ), 6) == 0)

# Test 2
p0 = c(6, 5)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 10.875/49.0625, 235.1875/49.0625 ), 6) == 0)

# Test 3
p0 = c(8, 3)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( -80/49.0625, 130.3125/49.0625 ), 6) == 0)

# Test 4
p0 = c(3, 8)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 147.1875/49.0625, 392.5/49.0625 ), 6) == 0)

# Test 5
p0 = c(2, 2)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 217.125/49.0625, 102.375/49.0625 ), 6) == 0)

cat("Finished Subsection [Reflections]: Slightly Perturbed Vertical Reflection - part 6\n");


###### Slightly Perturbed Horizontal Reflection ######

# Line Start = (0, 4), End = (5, 4.5)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(4, 4.5);
lines(x, y, lwd=2, col="red")
abline(v=x, col="green", lty=3)

# points(x, y)!
# Test 1
p0 = c(2, 6)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 59.5/25.25, 61.5 /25.25 ), 6) == 0)

# Test 2
p0 = c(1, 8)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 44.75/25.25, 7/25.25 ), 6) == 0)

# Test 3
p0 = c(3, 5)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 79.25/25.25, 91.25/25.25 ), 6) == 0)

# Test 4
p0 = c(4, 9)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 124/25.25, -2.75/25.25 ), 6) == 0)

# Test 5
p0 = c(4, 0)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 79/25.25, 220/25.25 ), 6) == 0)

# Test 6
p0 = c(0, 4)
p = testReflection(p0, x, y)
stopifnot(all(p == c(0,4)))

# Test 7
p0 = c(6, 2)
p = testReflection(p0, x, y)
testRoundValues(x, y, p)
stopifnot(round(p-c( 138.5/25.25, 180.5/25.25 ), 6) == 0)

cat("Finished Subsection [Reflections]: Slightly Perturbed Horizontal Reflection - part 7\n");

cat("Finished Section: Reflections \n")

##### Test: Shift point along line #####

###### Horizontal line ######

# Line Start = (0, 4), End = (5, 4)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(4, 4);
lines(x, y, lwd=2, col="red")

# Test 1 - origin of line
p0 = c(0, 4)
testShiftPoint(p0, x, y)

# Test 2 - over the line
p0 = c(2, 6)
testShiftPoint(p0, x, y, d=seq(1, 4, by=0.5), color="blue")

# Test 3 - under the line
p0 = c(2, 1)
testShiftPoint(p0, x, y)

# Test 4 - over the line
p0 = c(3, 7)
testShiftPoint(p0, x, y, d=seq(1,4, by=0.5), color="blue")

# Test 5 - shift behind
p0 = c(4, 4)
p = testShiftPoint(p0, x, y, d=-1)

cat("Finished Subsection [Shift point]: Horizontal line - part 1\n");

###### Origin: (0, y) ######

# Line Start = (0, 0), End = (5, 4)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(0, 4);
lines(x, y, lwd=2, col="red")

# Test 1 - origin of line
p0 = c(0, 0)
p = testShiftPoint(p0, x, y)

# Test 2 - over the line
p0 = c(2, 6)
p = testShiftPoint(p0, x, y, d=seq(1, 4, by=0.5), color = "blue")

# Test 3 - under the line
p0 = c(2, -1)
p = testShiftPoint(p0, x, y)

# Test 4 - over the line
p0 = c(3, 2)
p = testShiftPoint(p0, x, y, d=seq(1,4, by=0.5), color="blue")

# Test 5 - shift behind
p0 = c(4, 2)
p = testShiftPoint(p0, x, y, d = -1)

cat("Finished Subsection [Shift point]: Origin: (0, y) line - part 2\n");

###### Vertical line ######

# Line Start = (3, 1), End = (3, 8)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(3, 3); y = c(1, 8);
lines(x, y, lwd=2, col="red")

# Test 1 - origin of line
p0 = c(3, 8)
p = shiftPoint(p0, x, y, d=1)
points(p0[1], p0[2])
points(p, col="green")

# Test 2 - over the line
p0 = c(2, 2)
p = shiftPoint(p0, x, y, d=seq(1, 4, by=0.5))
points(p0[1], p0[2])
points(p, col="blue")

# Test 3 - under the line
p0 = c(4, 2)
p = shiftPoint(p0, x, y, d=1)
points(p0[1], p0[2])
points(p, col="green")

# Test 4 - over the line
p0 = c(0, 2)
p = shiftPoint(p0, x, y, d=seq(1,4, by=0.5))
points(p0[1], p0[2])
points(p, col="blue")

# Test 5 - shift behind
p0 = c(3, 2)
p = shiftPoint(p0, x, y, d=-1)
points(p0[1], p0[2])
points(p, col="green")

cat("Finished Subsection [Shift point]: Vertical line - part 3\n");

###### Ascending line ######
# Line Start = (0, 1), End = (10, 5) ##initial

plot.base()
x = c(0, 10); y = c(1, 5);
lines(x, y, lwd=2, col="red")

# Test 1 - origin of line
p0 = c(0, 1)
p = testShiftPoint(p0, x, y)

# Test 2 - over the line
p0 = c(0, 2)
p = testShiftPoint(p0, x, y, d=seq(1, 4, by=0.5), color = "blue")

# Test 3 - under the line
p0 = c(4, 2)
p = testShiftPoint(p0, x, y)

# Test 4 - over the line
p0 = c(0, 4)
p = testShiftPoint(p0, x, y, d=seq(1,4, by=0.5), color = "blue")

# Test 5 - shift behind
p0 = c(10, 4)
p = testShiftPoint(p0, x, y, d = -1)

cat("Finished Subsection [Shift point]: Ascending line - part 4\n");

###### Descending line ######

# Line Start = (0, 8), End = (8, 1)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 8); y = c(8, 1);
lines(x, y, lwd=2, col="red")

# Test 1 - origin of line
p0 = c(0, 8)
p = testShiftPoint(p0, x, y)

# Test 2 - over the line
p0 = c(6, 4)
p = testShiftPoint(p0, x, y, d=seq(1, 4, by=0.5), color = "blue")

# Test 3 - under the line
p0 = c(2, 4)
p = testShiftPoint(p0, x, y)

# Test 4 - over the line
p0 = c(4, 7)
p = testShiftPoint(p0, x, y, d=seq(1,4, by=0.5), color = "blue")

# Test 5 - shift behind
p0 = c(6, 2)
p = testShiftPoint(p0, x, y, d = -1)

cat("Finished Subsection [Shift point]: Descending line - part 5\n");

###### Slightly Perturbed Vertical Reflection ######

# Line Start = (3.25, 3), End = (1, 8)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(3.25, 3); y = c(1, 8);
lines(x, y, lwd=2, col="red")

# Test 1 - origin of line
p0 = c(3, 8)
p = testShiftPoint(p0, x, y)

# Test 2 - over the line
p0 = c(6, 2)
p = testShiftPoint(p0, x, y, d=seq(1, 4, by=0.5), color = "blue")

# Test 3 - under the line
p0 = c(5, 4)
p = testShiftPoint(p0, x, y)

# Test 4 - over the line
p0 = c(2, 3)
p = testShiftPoint(p0, x, y, d=seq(1,4, by=0.5), color = "blue")

# Test 5 - shift behind
p0 = c(3.5, 2)
p = testShiftPoint(p0, x, y, d = -1)

cat("Finished Subsection [Shift point]: Slightly Perturbed Vertical line - part 6\n");

###### Slightly Perturbed Horizontal Reflection ######

# Line Start = (0, 5), End = (4, 4.5)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(4, 4.5);
lines(x, y, lwd=2, col="red")

# Test 1 - origin of line
p0 = c(0, 4)
p = testShiftPoint(p0, x, y)

# Test 2 - over the line
p0 = c(2, 5)
p = testShiftPoint(p0, x, y, d=seq(1, 4, by=0.5), color = "blue")

# Test 3 - under the line
p0 = c(2, 3)
p = testShiftPoint(p0, x, y)

# Test 4 - over the line
p0 = c(3, 6)
p = testShiftPoint(p0, x, y, d=seq(1,4, by=0.5), color = "blue")

# Test 5 - shift behind
p0 = c(4, 4.25)
p = testShiftPoint(p0, x, y, d = -1)

cat("Finished Subsection [Shift point]: Slightly Perturbed Horizontal line - part 7\n");

cat("Finshed Section: Shift point along line \n")

##### Test: Shift Line #####

###### Horizontal line ######

# Line Start = (0, 4), End = (5, 4)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(4, 4);
lines(x, y, lwd=2, col="red")

# Test 1
testShiftLine(x, y)

cat("Finished Subsection [Shift line]: Horizontal line - part 1\n");

###### Origin: (0, y) ######

# Line Start = (0, 0), End = (5, 4)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(0, 4);
lines(x, y, lwd=2, col="red")

# Test 1
testShiftLine(x, y)

cat("Finished Subsection [Shift line]: Origin: (0, y) line - part 2\n");

###### Vertical line ######

# Line Start = (3, 1), End = (3, 8)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(3, 3); y = c(1, 8);
lines(x, y, lwd=2, col="red")

# Test 1
testShiftLine(x, y)

cat("Finished Subsection [Shift line]: Vertical line - part 3\n");

###### Ascending line ######
plot.base()
x = c(0, 10); y = c(1, 5);
lines(x, y, lwd=2, col="red")

# Test 1
testShiftLine(x, y)

# Line Start = (0, 8), End = (8, 1)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 8); y = c(8, 1);
lines(x, y, lwd=2, col="red")

# Test 2
testShiftLine(x, y)

cat("Finished Subsection [Shift line]: Ascending line - part 4\n");

###### Slightly Perturbed Vertical line ######

# Line Start = (3.25, 1), End = (3, 8)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(3.25, 3); y = c(1, 8);
lines(x, y, lwd=2, col="red")

# Test 1
testShiftLine(x, y)

cat("Finished Subsection [Shift line]: Slightly Perturbed Vertical line - part 5\n");

###### Slightly Perturbed Horizontal line ######

# Line Start = (0, 4), End = (5, 4.5)
plot.base(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2))
x = c(0, 5); y = c(4, 4.5);
lines(x, y, lwd=2, col="red")

# Test 1
testShiftLine(x, y)

cat("Finished Subsection [Shift line]: Slightly Perturbed Horizontal line - part 6\n");

cat("Finished Section: Shift Line \n")


###############
### Circles ###

### Circle Intersections
center1 = c(2,3)
alpha = pi/5
r = 3/2
# Center 2 is ON circle 1;
mid2 = center1 + r * c(cos(alpha), sin(alpha));
# p = solve.circle.intersection2(mid2, center1, 0.5, r)
xy = rbind(center1, mid2);
p = solve.circle.intersection(xy[,1], xy[,2], c(r, 0.5))
plot.base()
shape::plotellipse(r, r, mid=center1)
shape::plotellipse(0.5, 0.5, mid=mid2, lcol="red")
points(p$x, p$y, col="green", lwd=2)
id = c(2,1)
stopifnot(round(p$x - c(2.856325, 3.435889)[id], 5) == 0)
stopifnot(round(p$y - c(4.231547, 3.433845)[id], 5) == 0)

