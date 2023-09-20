


### Cholesterol

#' @export
cholesterol = function(origin = c(0,0), dPi = 0.125) {
	lst = cholesterol.bb(origin = origin);
	lst$Me1 = ligand(lst, 1, 1, slope = Inf);
	lst$Me2 = ligand(lst, 3, 1, slope = Inf);
	lst$Pi  = as.piBond(lst[2], which = c(4,5), d = - dPi)[[1]];
	# OH
	lstOH = ligand(lst, 1, i = 4, slope = 1/sqrt(3), d = -0.75)[[1]];
	txtOH = list(x = lstOH$x[2], y = lstOH$y[2], labels = "HO", adj = c(1.2, 0.8));
	class(txtOH) = c("text", "list");
	lstOH = list(OHb = lstOH, OHl = txtOH);
	class(lstOH) = c("bioshape", "list");
	lst$OH = lstOH;
	# TODO: 17-side chain;
	return(lst);
}
#' @export
cholesterol.bb = function(origin = c(0,0)) {
	lst = parseCycles("6|6\\6|5");
	if(any(origin != 0)) {
		lst = shift.listNR(lst, by = origin);
	}
	return(lst);
}

#' @export
example.cholesterol = function(origin = c(0,0)) {
	par.old = par(mfrow = c(2,1))
	
	### Cholesterol backbone
	plot.base(ylim = c(-2,4) + origin[2]);
	lines(cholesterol.bb(origin = origin));
	
	### Cholesterol
	plot.base(ylim = c(-2,4) + origin[2]);
	lines(cholesterol(origin = origin));
	
	par(par.old);
}


######################

### Polycycles

#' @export
test.polycycle = function(n = c(14, 17, 11, 19), ngon = 7,
		Ht = c("O", "S", NA, ""), col = c("red", "#FFE832", "grey", "grey"), R = 4) {
	if(length(R) == 1) R = rep(R, 4);
	if(length(ngon) == 1) ngon = rep(ngon, 4);
	#
	par.old = par(mfrow = c(2,2));
	for(i in seq(4)) {
		lim = R[i] + 1; lim = c(- lim, lim);
		plot.base(xlim = lim, ylim = lim, axt = NULL, asp = 1);
		if(is.na(Ht[i])) {
			lines(polycycle.polyene(n[i], ngon[i]));
		} else if(nchar(Ht[i]) == 0) {
			lines(polycycle.polyEE(n[i], ngon[i]));
		} else {
			lines(polycycle.polyHt(n[i], ngon[i], Ht = Ht[i], col = col[i]));
		}
		lines(annotate.polycycle(n[i]))
	}
	par(par.old);
	invisible();
}


###########################

### Test: Newman Projection

#' @export
test.proj.newman = function(phi.add = c(0,0), ligand = "F,H,H|OH,\\CH[3],H") {
	plot.base()
	# Row 1:
	lines(proj.newman(ligand, center = c(1,6), phi = c(0.7, 1) + phi.add))
	lines(proj.newman(ligand, center = c(6,6), phi = c(1.6, 1.3) + phi.add))
	# Row 2:
	lines(proj.newman(ligand, center = c(1,1), phi = c(0.5, 0.8) + phi.add))
	lines(proj.newman(ligand, center = c(6,1) + phi.add))
}
