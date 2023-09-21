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


####################

### Helper Functions

### Rep List
#' @export
rep.all.list = function(x, len) {
	if(len == 1) return(x);
	ll = length(x);
	for(id in seq(ll)) {
		tmp = x[[id]];
		if(length(tmp) == 1) x[[id]] = rep(tmp, len);
	}
	return(x);
}

#' @export
merge.first.list = function(x, y) {
	d = setdiff(y, x);
	c(x, d);
}


# xy    = list(x, y);
# new   = list(x, y) of matrixes with new values;
# which = ids (in x & y) to replace;
#' @export
replace.byVal = function(xy, new, which) {
	len = length(xy);
	if(len == 0 || length(which) == 0) return(xy);
	for(id in seq(len)) {
		tmp = xy[[id]];
		tmp$x[which] = new$x[, id];
		tmp$y[which] = new$y[, id];
		xy[[id]] = tmp;
	}
	return(xy);
}


### Rbind with R-Shifted Matrix
#' @export
rbind.shiftR = function(x) {
	len = ncol(x);
	tmp = cbind(x[, len], x[, -len]);
	rbind(x, tmp);
}


### Drop column
#' @export
drop.col = function(x, name) {
  id = match(name, names(x));
  isNA = is.na(id);
  if(any(isNA))
    warning("Names not found: ", paste0(name[isNA], collapse=", "));
  id = id[ ! isNA];
  if(length(id) == 0) return(x);
  return(x[, -id]);
}


##################

### Rolling Apply:
#' @export
apply.roll = function(x, FUN, ...) {
	len = length(x);
	id  = seq(len);
	idc = c(id[-1], 1); # cyclic shift left: Next element;
	rr  = lapply(id, function(id) FUN(x[[id]], x[[idc[id]]], ...));
	do.call(rbind, rr);
}


### Replace with Midline
# xy = list of elements of type list(x, y);
#' @export
as.mean.xy = function(xy, which1, which2) {
	mean.xy2 = function(e1, e2, which1, which2) {
		x = (e1$x[which1] + e2$x[which2]) / 2;
		y = (e1$y[which1] + e2$y[which2]) / 2;
		data.frame(x = x, y = y);
	}
	mxy = apply.roll(xy, mean.xy2, which1=which1, which2=which2);
	x = matrix(mxy$x, nrow=2);
	y = matrix(mxy$y, nrow=2);
	mxy = list(x = rbind.shiftR(x), y = rbind.shiftR(y));
	replace.byVal(xy, mxy, which = c(which1, which2));
}

##################

### Formulas

### Negate Formula
#' @export
formula.neg = function(x) {
	tmp = expression(- (0))[[1]];
	tmp[[2]][[2]] = x[[2]];
	x[[2]] = tmp;
	return(x)
}

