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


#### Tests ####

### Simple Circles
# - plot radially n circles;
n = 10
R = 8
r = 1
phi = pi / n; # add some rotation
# Ex 1a:
xy = pointsCircle(n, R, phi=phi);
test.FilledCircle(xy, r=r, lim = R+1, lcol=1);
# Ex 1b:
test.FilledCircle(xy, r=seq(r,2, length.out=n), lim = R+1, lcol=1);


#####################

#####################
### Circle Chains ###

### Closed Circles
# - n & r known, R unknown;
n = 15
r = 1
phi = pi / n; # add some rotation
xy = circles.OnCircle(n, r, phi=phi);
test.FilledCircle(xy);


### Radius of Big Circle: known
# - small circles: on the big circle;
# - n & R known, r unknown;
n = 15
R = 7
xy = circles.OnFixedCircle(n, r=R, phi=phi);
test.FilledCircle(xy, R=R);


#### Inside Circle ####

# Tangent & Inside given Circle
# Outer Circle: known

### Example 1: Simple
n = 19
R = 15
phi = pi / n; # add some rotation
xy = circles.InFixedCircle(n, R, phi=phi);
test.FilledCircle(xy, R=R)


### Example 2:
# - Circle 1: R unknown
n = 19
r = 1
phi = pi / n; # add some rotation
xy = circles.OnCircle(n, r, phi=phi);
test.FilledCircle(xy)

# - Circle 2: based on R of Circle 1;
# Outer Circle: known
# Inner Circle: unknown
R = attr(xy, "R") - r;
xy = circles.InFixedCircle(n, r=R, phi=phi);
test.FilledCircle(xy, add=TRUE, line=FALSE);


#### Outside Circle ####

# Tangent & Outside given Circle
n = 13
R = 6
phi = pi / n; # add some rotation
xy = circles.OutsideFixedCircle(n, R, phi=phi);
test.FilledCircle(xy, R=R);


#### Shifted Center ####

# 2 Circles
n = 23
R = 6
phi  = pi / n;
mid1 = c(-R, 0); mid2 = mid1 + c(2*R, 0);
xy1 = circles.InFixedCircle(n, r=R, center=mid1);
xy2 = circles.InFixedCircle(n, r=R, center=mid2, phi=phi);
test.FilledCircle(xy1, R=R, lim = 2*R + 1);
test.FilledCircle(xy2, R=R, add=TRUE);


#####################

#### Cell Shapes ####

### Brush-Border Cell
p1 = c(0,0)
cell = cell.BrushBorder(p1, w=5, h=8)
plot.base()
lines(cell)


### Slanted Cell
p1 = c(5,8); slope = -1;
cell = cell.BrushBorder(p1, w=5, h=8, slope=slope)
plot.base()
lines(cell)


### Smooth Muscles / Connective Tissue
plot.base()
lst = cell.SmoothMuscle(c(1,8), c(3, 7))
lines(lst, lwd=2)

lst = cell.SmoothMuscle(c(2,8), c(1, 5), r=0.6)
lines(lst, lwd=2)
abline(v = c(1,8), col = "green")
abline(h = c(1,7), col = "green")


##################

#### Liposome ####

plot.base(xlim=c(-10,10), ylim=c(-10,10))
lst = liposomes(c(30, 17), r=0.5, phi=c(0, pi/34), d=0.2)
lines(lst)


### Helix / DNA

### Ex 1a:
p1 = c(1,1); p2 = c(8,3);
lst1 = helix(p1, p2)
lst2 = helix(p1, p2, phi=-pi/2)
plot.base()
lines(lst1)
lines(lst2)
# visual aids
abline(v = p1[1], col="green")
abline(v = p2[1], col="green")

### Ex 1b: Shifted version of helices
dy = c(0, 3); n = 2.5;
lst1 = helix(p1 + dy, p2 + dy, n=n)
lst2 = helix(p1 + dy, p2 + dy, n=n, phi=-pi/2)
# plot.base() # add to previous image;
lines(lst1)
lines(lst2)

### Ex 1c: Reversed direction
# different p2!
dy = c(0, 5); p2 = c(7,9)
lst1 = helix(p1 + dy, p2)
lst2 = helix(p2, p1 + dy, phi=-pi/2)
# plot.base() # add to previous image;
lines(lst1, col="blue", lwd=2)
lines(lst2, col="red", lwd=2)


### Ex 2: Vertical Object
p1 = c(1,1); p2 = c(p1[1], 8);
lst1 = helix(p1, p2)
lst2 = helix(p1, p2, phi=-pi/2)
plot.base()
lines(lst1)
lines(lst2)
# visual aids:
abline(h = p1[2], col="green")
abline(h = p2[2], col="green")

# reversed direction:
dx = c(4, 0)
lst1 = helix(p1 + dx, p2 + dx)
lst2 = helix(p2 + dx, p1 + dx, phi=-pi/2)
lines(lst1, col="blue")
lines(lst2, col="red")



#### Spirals / Coils ####

### Ex 1:
p1 = c(1,1); p2 = c(2,8); dx = c(2.5,0);
lst1 = spirals(p1, p2)
lst2 = spirals(p1 + dx, p2 + dx)
plot.base()
lines(lst1)
lines(lst2)

### Ex 2:
p1 = c(1,1); p2 = c(2,8); dx = c(3,0);
lst1 = spirals(p1, p2)
lst2 = spirals(p2 + dx, p1 + dx)
plot.base()
lines(lst1)
lines(lst2)

### Ex 3: Vertical
p1 = c(1,1); p2 = c(p1[1],8); dx = c(2.5,0);
lst1 = spirals(p1, p2)
lst2 = spirals(p2 + dx, p1 + dx)
plot.base()
lines(lst1)
lines(lst2)

