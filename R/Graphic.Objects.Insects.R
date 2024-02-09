

### Fly
# see eg. Image 1 in Vignette:
# https://cran.r-project.org/web/packages/metaRange/vignettes/A01_intro.html
# bhsift = overlap Body with Head;
#' @export
fly = function(x, y, r = c(1.5, 2, 3), rw = c(1, 3), lwd = 2,
		col = 1, fill = c("white", "#D0D0A0", "white"),
		alpha = 75, bshift = 0.25) {
	if(length(lwd) == 1) lwd = rep(lwd, 3);
	cH = c(x[1], y[1]);
	cB = c(x[1], y[1] - r[1] * (1 - bshift) - r[3]); # TODO
	H  = as.circle(list(center = cH, r = r[1], lwd = lwd[1], col=col, fill = fill[1]));
	B  = ellipse(cB, r = r[c(2,3)], lwd = lwd[2], col=col, fill = fill[2]);
	### Wings
	phi = alpha * pi / 360; # half of angle between Wings;
	dx  = rw[2] * cos(phi) + r[1] / 3; # hardcoded: 1/3;
	dy  = rw[2] * sin(phi);
	yW  = y[1] - r[1] - r[3] / 2 - dy; # hardcoded: 1/2;
	lwd = lwd[3];
	W1 = ellipse(x[1] + dx, yW, r = rw, theta = phi,
		lwd=lwd, col=col, fill = fill[3]);
	W2 = ellipse(x[1] - dx, yW, r = rw, theta = - phi,
		lwd=lwd, col=col, fill = fill[3]);
	### Fly
	lst = list(Head = H, Body = B, Wings = as.bioshape(list(W1, W2)));
	return(as.bioshape(lst));
}

