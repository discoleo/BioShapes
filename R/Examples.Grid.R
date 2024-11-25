#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes


### Hexagonal Grid: Rings
# phiX = Rotation of last set of hexagonal grids;
#' @export
example.grid.hexa = function(n = c(4,5,5,5), center = c(4,4),
		phiX = c(0,-1,1,-2,2) * pi / 12,
		col = c("red", "#32D064B0", "#3280D4B0", "#E644A0B0", "#E6A460B0")) {
	par.old = par(mfrow = c(2,2));
	# Ex 1:
	r = 1/2; phi = 0;
	hg = grid.hexa(n[1], r=r, center=center, col=col, phi=phi)
	plot.base()
	lines(hg)
	# Ex 2:
	r = 1/2; phi = pi/10;
	hg = grid.hexa(n[2], r=r, center=center, col=col, phi=phi)
	plot.base(axt = NULL)
	lines(hg)
	# Ex 3:
	r = 1/2; phi = pi/10;
	plot.base(axt = NULL)
	hg = grid.hexa(n[3], r=r, center=center, col=col, phi=phi)
	lines(hg)
	hg = grid.hexa(n[3], r=r, center=center + c(0.25,0), col=col, phi=phi)
	lines(hg)
	# Ex 4:
	r = 1/2;
	plot.base(axt = NULL)
	for(phi in phiX) {
		hg = grid.hexa(n[4], r=r, center=center, col=col, phi=phi);
		lines(hg);
	}
	#
	par(par.old);
	invisible();
}