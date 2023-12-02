#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#   GitHub: https://github.com/Adi131313/BioShapes
#
# 2. Bachelor Thesis: Darian Voda (2021-2022)


### Curves

#' @export
example.curves = function(R = 5, nr = 20, axt = c(1, 2)) {
  par.old = par(mfrow = c(2,2));

  ### Circular Helix
  lim = c(-R, R) + c(-1,1);
  plot.base(xlim = lim, ylim = lim)
  xy = helix.rad(R=R, n=nr, N=512)
  lines(xy)
  xy = helix.rad(R=R, n=nr, phi=pi/2, N=512)
  lines(xy, col="red")


  ### Flower 1:
  lim = c(-R, R) + c(-2,2);
  plot.base(xlim = lim, ylim = lim)
  xy = helix.rad(R=R, n=nr-9, phi=pi)
  lines(xy)
  xy = helix.rad(R = R - 1, n=nr-9, phi=pi)
  lines(xy, col="red")


  ### Flower 2:
  lim = c(-R, R)*2 + c(-1,1);
  plot.base(xlim = lim, ylim = lim)
  xy = helix.rad(R=R, r=R, n=round(nr/2), N=256)
  lines(xy, col="red")


  ### Other
  lim = c(-R, R) + c(-2,2);
  plot.base(xlim = lim, ylim = lim)
  abline(h=0, col="green", lty=2)
  for(r in c(3,4,6,9)) {
    xy = helix.rad(R = R - 2/r, n=nr-10, phi=pi, r = 2/r);
    lines(xy, col="red");
  }
  par(par.old);
  invisible();
}

#########################

#########################
### DNA: Double Helix ###

### Basic Example
#' @export
example.dna = function(n = 3, n.lin = 6, phi = c(0, pi),
		lwd = 2, lwd.lines = lwd, col = c("red", "green"), axt = NULL){
	h2 = dna.new(c(1,8), c(1,4), phi = phi,
		n=n, n.lin=n.lin, lwd=lwd, lwd.lines=lwd.lines, col=col);
	plot.base(axt=axt);
	lines(h2);
	invisible(h2);
}

### Variants: dphi
# dphi = shift of phi;
# both.phi = shift phi for both strands;
#' @export
example.dna.tests = function(n = 3, phi = c(pi, -pi/2), A = 1,
		x0 = c(0, 8), y0 = c(0, 0), dy = 3,
		dphi = c(0, pi/6, pi/3, pi/2), both.phi = TRUE,
		col = c("red", "green"),
		xlim = NULL, ylim = NULL, axt = c(1,2), verbose = FALSE) {
	if(length(n) > 1) stop("Invalid number of helix turns!");
	LEN = length(dphi);
	if(length(dy) == 1) {
		dy = seq(0, LEN - 1) * dy;
	}
	if(is.null(xlim)) {
		sgn  = if(x0[1] < x0[2]) c(-1, 1) else c(1, -1);
		xlim = x0 + 2*sgn;
	}
	if(is.null(ylim)) ylim = c(-2, 10);
	plot.base(xlim=xlim, ylim=ylim, axt=axt);
	for(i in seq(LEN)) {
		phi_i = phi + if(both.phi) dphi[i] else c(0, dphi[i]);
		if(verbose) cat("Phi = ", phi_i / pi, "\n");
		lines(dna.new(x0, y0 + dy[i], A=A, col=col, n=n, phi = phi_i));
	}
}

# more Tests:
#' @export
test.dna.nb = function(phi = pi, n = c(3, 3.4, 3, 1.7), y = list(c(0, -2))) {
	if(length(phi) == 1) phi = c(0, phi);
	par.old = par(mfrow = c(2,2))
	###
	example.dna.tests(n = n[1], phi=phi)
	
	###
	example.dna.tests(n = n[2], phi=phi)

	###
	example.dna.tests(n = n[3], phi=phi, y0 = y[[1]])
	
	###
	example.dna.tests(n = n[4], phi=phi, x0 = c(0, 6), y0 = y[[1]])
	# Note: (x0, y0) are NOT exact;
	abline(v=6, col="blue")
	
	par(par.old)
}

### Test: Bug in Function helix
#' @export
test.helix.directions = function() {
	plot.base();
	p1 = c(3,1); p2 = c(1,5)
	xy = helix(p1, p2); lines(xy);
	points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red")
	#
	p1 = c(4,6); p2 = c(1, 5)
	xy = helix(p1, p2); lines(xy);
	points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red")
	#
	p1 = c(4,6); p2 = c(6,1)
	xy = helix(p1, p2); lines(xy);
	points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red")
	#
	p2 = c(4,6); p1 = c(1, 5)
	xy = helix(p1, p2); lines(xy);
	points(c(p1[1], p2[1]), c(p1[2], p2[2]), col="red")
}


##################
##################

### Connected-Arcs

#' @export
test.curve.c2pi = function(x = c(1, 6), y = c(3, 6), r = c(2, 1.5),
		col = c("black", "blue")) {
	cc = solve.curve.c2pi(x, y, r);
	plot.base()
	if(is.na(cc$C1[[1]])) return(cc); # NO solution;
	lines(as.circle(list(center = cc$C1, r = r[1], col = col[[1]])))
	lines(as.circle(list(center = cc$C2, r = r[2], col = col[[1]])))
	points(x, y, col = "red")
	#
	cc = solve.curve.c2pi(x, y, r, rev = TRUE);
	lines(as.circle(list(center = cc$C1, r = r[1], col = col[[2]])))
	lines(as.circle(list(center = cc$C2, r = r[2], col = col[[2]])))
}