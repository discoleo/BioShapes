#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes


### Helper Functions: Operations on Circles


### Intersection: 2 Circles
#' @export
solve.circle.intersection2 = function(center1, center2, r1, r2, digits=4, debug=TRUE) {
  xp = center1[[1]]; yp = center1[[2]];
  xc = center2[[1]]; yc = center2[[2]];
  d  = r1; r = r2;
  xs = xp + xc;
  xd = xp - xc; xd2 = xs*xd;
  yd = yp - yc;
  # P2:
  a  = 4*(xd^2 + yd^2);
  b1 = - 4*(xs^3 - 4*xp*xc*xs + yd^2*xs + (r^2 - d^2)*xd);
  b0 = d^4 + r^4 + xc^4 + xp^4 + yp^4 + yc^4 - 4*yp^3*yc - 4*yp*yc^3 + 6*yp^2*yc^2 +
    + 2*xp^2*yp^2 + 2*xp^2*yc^2 + 2*xc^2*yp^2 + 2*xc^2*yc^2 - 2*xp^2*xc^2 +
    + 2*(r^2 - d^2)*xd2 - 2*(r^2 + d^2)*yd^2 +
    - 4*yp*yc*xp^2 - 4*yp*yc*xc^2 - 2*d^2*r^2;
  # Det:
  D = b1^2 - 4*a*b0;
  if(debug) { cat("Coeffs: "); print(c(a, b1, b0, D)); }
  if(round(D, digits) == 0) D = 0;
  if(D < 0) return(NULL);
  D = sqrt(D);
  # Sol:
  x = ( - b1 + c(-D, D)) / (2*a);
  y =  yp^2 - yc^2 + xd2 - 2*xd*x + r^2 - d^2;
  div = 2*yd; # TODO: yp == yc;
  y = y / div;
  sol = data.frame(x=x, y=y);
  return(sol);
}

# simplified approach
#' @export
solve.circle.intersection = function(x, y, r) {
	dL = dist.xy(x, y);
	dr = (r[1] + r[2] - dL) / 2;
	if(dr < 0) return(NA);
	ss = (r[1] + r[2] + dL) / 2;
	area = sqrt(ss * dr * (ss - r[1]) * (ss - r[2]));
	slope = slope(x, y);
	h  = 2 * area / dL;
	d1 = sqrt(r[1]^2 - h^2);
	t1 = d1 / dL;
	xc = (1 - t1)*x[1] + t1*x[2];
	yc = (1 - t1)*y[1] + t1*y[2];
	xy = shift.ortho(xc, yc, d = c(-h, h), slope=slope);
	return(xy);
}

