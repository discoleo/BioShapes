


###############
### Circles ###

test.isNaN = function(x) {
	if(is.nan(x) || is.na(x)) {
		return(FALSE);
	}
	return(TRUE);
}

### Constructors

### Through 3 Points

test.circle.p3 = function(x, y, msg="") {
	cc = circle.p3(x, y)
	test.r = function(x) {
		test.isNaN(x$r);
	}
	if( ! test.r(cc)) {
		cat(paste0("Failed in ", msg, ": x = ",
			x[1], ", ", x[2], "; y = ", y[1], ", ", y[2], ";\n"));
		stop();
	}
	#
	plot.base()
	lines(as.bioshape(list(cc)))
	polygon(x, y, border="blue")
	points(x, y, col="red")
}

###
par.old = par(mfrow = c(3,3))


# VH
x = c(1,1,5); y = c(1,4,1)
test.circle.p3(x, y)

# VH
x = c(1,6,1); y = c(1,1,4)
test.circle.p3(x, y)


# V-Narrow
x = c(1,1,5); y = c(1,4,3)
test.circle.p3(x, y)

# V-Obtuse
x = c(1,1,5); y = c(1,4,0)
test.circle.p3(x, y)

# V-Obtuse
x = c(7,1,7); y = c(4,2,8)
test.circle.p3(x, y)

# H-Narrow
x = c(1,5,4); y = c(1,1,8)
test.circle.p3(x, y)

# H-Obtuse
x = c(1,8,6); y = c(1,1,0)
test.circle.p3(x, y)

# H-Obtuse
x = c(7,1,3); y = c(8,2,8)
test.circle.p3(x, y)

# G-Obtuse
x = c(1,8,6); y = c(1,2,0)
test.circle.p3(x, y)

###
par(par.old)


##############################
##############################

### Tangent & through 2 Points

test.circle.p2s()
test.circle.p2s(slope = -2, add=T, col = "#A032F0", col.line = "#88A0D2")


#################

### Arcs / Curves

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


