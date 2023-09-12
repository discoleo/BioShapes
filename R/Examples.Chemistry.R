


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