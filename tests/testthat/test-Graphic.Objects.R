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

