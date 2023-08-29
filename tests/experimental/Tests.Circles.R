



### Circles

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


