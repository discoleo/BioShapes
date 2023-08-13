###################
#
# Title: BioShapes
#
# Maintainer: L. Mada
#
#
# Continues the work of:
# - Adrian Cotoc (BSc 2022-2023)
# - Darian Voda (BSc 2021-2022)
#
# GitHub: https://github.com/discoleo/BioShapes


### Arrows: Other Types

# - Pins/Tags;

#######################

### Pin/Tag
# (x, y) = from (pin-point, center);
# theta  = angle of pin;
#' @export
pin.center = function(x, y, theta = pi/3, lwd=1, lwd.circle = lwd,
		col = 1, col.circle = col, fill = NULL, debug = FALSE) {
	th = theta / 2;
	d = dist.xy(x, y);
	r = d * sin(th);
	# Circle Arc
	slope = slope(x, y);
	phi = atan2(y[2] - y[1], x[2] - x[1]);
	phi = - pi/2 + (phi - th);
	phi = c(phi, phi - (pi - theta));
	if(debug) print(phi/pi);
	lst = list(r = r, center = c(x[2], y[2]), phi = phi,
		lwd = lwd.circle, col = col.circle);
	if( ! is.null(fill)) lst$fill = fill;
	class(lst) = c("circle.arc", "list");
	# Triangle:
	tt = cos(th)^2;
	xx = tt*x[2] + (1-tt)*x[1];
	yy = tt*y[2] + (1-tt)*y[1];
	hT = d * sin(theta) / 2;
	pp = shift.ortho(c(xx, yy), slope = slope(x, y), d = c(-hT, hT));
	xx = c(pp[1,1], x[1], pp[2,1]);
	yy = c(pp[1,2], y[1], pp[2,2]);
	lstT = list(x = xx, y = yy, lwd=lwd, col=col);
	if( ! is.null(fill)) {
		lstP = lstT;
		lstP$col = NA; # remove border;
		lstP$fill = fill;
		class(lstP) = c("polygon", "list");
		lstT = list(Poly = lstP, Triangle = lstT);
		lstT = as.bioshape(lstT);
	}
	lst0 = list(Triangle = lstT, Arc = lst);
	return(as.bioshape(lst0));
}

