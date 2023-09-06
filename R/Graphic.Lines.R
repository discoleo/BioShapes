

### Lines & Connectors


### Line with Stair
#' @export
as.lines.stairs = function(x, y, slope, t = 0.5) {
	d = dist.xy(x, y);
	s = sqrt(slope^2 + 1);
	h = (slope * (x[2] - x[1]) - (y[2] - y[1])) / s;
	D = sqrt(d^2 - h^2);
	qd = which.quadrant(x, y);
	sg = if(qd == 1 || qd == 4) 1 else -1;
	D  = sg * D;
	p1 = shift.point.default(c(x[1], y[1]), slope=slope, d = D);
	p2 = shift.point.default(c(x[2], y[2]), slope=slope, d = - D);
	# Breakpoint:
	mid1.x = (1-t)*x[1] + t*p1[1];
	mid1.y = (1-t)*y[1] + t*p1[2];
	mid2.x = (1-t)*p2[1] + t*x[2];
	mid2.y = (1-t)*p2[2] + t*y[2];
	x = c(x[1], mid1.x, mid2.x, x[2]);
	y = c(y[1], mid1.y, mid2.y, y[2]);
	lst = list(x=x, y=y);
	return(lst);
}

### Stair: Equal Size
#' @export
as.lines.stairsEq = function(x, y, up = TRUE) {
	d = dist.xy(x, y);
	h = d / sqrt(5);
	#
	qd = which.quadrant(x, y);
	sg = if(qd == 1 || qd == 4) 1 else -1;
	D  = 2 * sg * h;
	s0 = slope(x, y);
	slope = if(up) (2*s0 - 1) / ( 2 + s0)
		else (2*s0 + 1) / ( 2 - s0);
	p1 = shift.point.default(c(x[1], y[1]), slope=slope, d = D);
	p2 = shift.point.default(c(x[2], y[2]), slope=slope, d = - D);
	x1 = (x[1] + p1[1])/2; y1 = (y[1] + p1[2])/2;
	x2 = (x[2] + p2[1])/2; y2 = (y[2] + p2[2])/2;
	#
	x = c(x[1], x1, x2, x[2]);
	y = c(y[1], y1, y2, y[2]);
	lst = list(x=x, y=y);
	return(lst);
}

### Connect by Stair
#' @export
connect.lines.stairs = function(x, y, slope, t = 0.5) {
	lst = as.lines.stairs(x, y, slope=slope, t=t);
	lst = list(L = lst);
	return(as.bioshape(lst));
}

### Connect: Stair with Equal Size
# up = upward stair;
#' @export
connect.lines.stairsEq = function(x, y, up = TRUE) {
	lst = as.lines.stairsEq(x, y, up=up);
	lst = list(L = lst);
	return(as.bioshape(lst));
}

