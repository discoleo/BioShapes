


###############
### Circles ###

### Constructors

### Through 3 Points

test.circle.p3Tests()


##############################
##############################

### Tangent & through 2 Points

test.circle.p2s()
test.circle.p2s(slope = -2, add=T, col = "#A032F0", col.line = "#88A0D2")


#################

### Arcs / Curves

# 2xC passing through 2 points
R = c(2, 1.5); # R = c(1.5, 2);
x = c(1, 6); y = c(3, 6);
test.curve.c2pi(x, y, r = R)


#####################
#####################

### Uniform Points

# Spiral Roulette
test.circle.uniform.text(c(100,100), phi=c(0, - pi))

# Tumor Mass: 2 cell types
test.tumor.mass2()


################

################
### Ellipses ###

### Tangents at Ellipse
test.ellipse.tan(c(2,4.5), r=c(1, 3), phi = -pi/5, center = c(3,5))

test.ellipse.tan(c(2,4.5), r=c(1, 3), phi = pi/3, center = c(3,5))

# Intersections with Parallel Lines
test.ellipse.intersect()

# Construct Ellipses
test.ellipse.bySlope()


### Cylindeers

test.cylinder()

test.cylinder(lty.back = NULL)


#####################

### Circle Chains ###

n = 15
plot.base()
lines(circles.TanToChainShape(n, 4, center = c(4,4), phi = pi/n))


