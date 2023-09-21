



### Circles

test.r = function(x) {
	test.isNaN(x$r);
}
test.isNaN = function(x) {
	if(is.nan(x) || is.na(x)) {
		return(FALSE);
	}
	return(TRUE);
}

### Through 3 Points

test.circle.p3 = function(x, y, msg="") {
	cc = circle.p3(x, y)
	if( ! test.r(cc)) {
		cat(paste0("Failed in ", msg, ": x = ",
			x[1], x[2], "; y = ", y[1], y[2], ";\n"));
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
test.circle.p2s = function(p1 = c(1,8), p2=c(5,5), slope = 1/12, dx=c(-5, 5),
		col = 1, col.line = "#32D048", add = FALSE) {
	tmp = circle.p2s(p1, p2, slope = slope);
	tmp = list(tmp, col = col);
	#
	if( ! add) plot.base();
	lines(as.bioshape(tmp));
	for(dxi in dx)
		lines(c(p1[1], p1[1] + dxi), c(p1[2], dxi*slope + p1[2]), col=col.line)
	points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red");
	return(invisible(tmp));
}

test.circle.p2s()
test.circle.p2s(slope = -2, add=T)


#####################
#####################

### Hashed Circle

test.circle.hash = function(n, phi, r = 1, scale = 1, center = c(4,4), ...) {
	tmp = circle.hash(n, phi=phi, r=r, center=center, scale=scale, ...)
	plot.base(asp = 1/scale)
	lines(tmp)
	abline(v = center[1] + c(-r,r), col = "green")
}


###
par.old = par(mfrow = c(2,2))

### Base:
test.circle.hash(11, phi = pi/5, r = 2, center = c(4,4), lty = 2, col = "red")

### Scale
n = 13; r = 2; scale = 2;
test.circle.hash(n, phi = pi/3, r=r, scale=scale)
#
test.circle.hash(n, phi = pi/2 + 0.2, r=r, scale=scale)
#
phi = pi - 0.2
test.circle.hash(n, phi = phi, r=r, scale=scale)
lines(circle.hash(11, phi = phi + pi/2, r=r, scale=scale, center = c(4,4), col = "red"))

par(par.old)


################

################
### Ellipses ###

### Tangents at Ellipse
test.ellipse.tan(c(2,4.5), r=c(1, 3), phi = -pi/5, center = c(3,5))

test.ellipse.tan(c(2,4.5), r=c(1, 3), phi = pi/3, center = c(3,5))


### Cylindeers

test.cylinder()

test.cylinder(lty.back = NULL)

