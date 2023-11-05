#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. BSc Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#   [old] GitHub: https://github.com/Adi131313/BioShapes
#
# 2. BSc Thesis: Darian Voda (2021-2022)


### Arrows: Helper Functions

as.arrow = function(x) {
	if( ! inherits(x, "arrow")) {
		class(x) = c("arrow", class(x));
	}
	invisible(x);
}


### Arrow Tail:
#' @export
arrowTail = function(x, y, d.lines, lwd=1, slope=NULL, scale = 1) {
  if(is.null(slope)) slope = slope(x, y);
  if(any(d.lines != 0)) {
    arrTail = shift.ortho(x, y, d = d.lines, slope=slope, scale=scale);
  } else {
    arrTail = list(x=x, y=y);
  }
  arrTail = list(xy = arrTail, lwd=lwd);
  return(arrTail)
}

# Tail = simple lines;
# xy  = data.frame with the Tail;
# xyH = data.frame with the Head;
# type:
# - Real = true intersection;
# - In   = xy will be intersected by extension of xyH;
# - Out  = extension of xy will intersect xyH;
intersect.arrow = function(xy, xyH, type = 1) {
	if(is.character(type)) {
		type = pmatch("Real", "In", "Out");
		if(is.na(type)) stop("Invalid intersection type!");
	}
	idT = unique(xy$id);
	if(length(idT) == 0) return(xy);
	isList = TRUE;
	if(inherits(xyH, "data.frame")) {
		isList = FALSE;
		idH = unique(xyH$id);
	} else {
		len = length(xyH$x);
		idH = if(len > 0) seq(len - 1) else numeric(0);
	}
	if(length(idH) == 0) return(xy);
	isIntersect = if(type == 3) {
		function(tt) return(tt >= 0 && tt <= 1);
	} else function(tt) { return(FALSE); };
	xyT = lapply(idT, function(id) {
		nR = which(xy$id == id);
		x  = xy$x[nR]; y = xy$y[nR];
		for(idZ in idH) {
			if(isList) {
				idZ = c(idZ, idZ + 1);
				xHi = xyH$x[idZ]; yHi = xyH$y[idZ];
			} else {
				isZ = xyH$id == idZ;
				xHi = xyH$x[isZ]; yHi = xyH$y[isZ];
			}
			xyI = intersect.lines(x, y, xHi, yHi);
			if(is.intersect.lines(xyI) || isIntersect(xyI$t2[1])) {
				x[2] = xyI$x[1];
				y[2] = xyI$y[1];
				xy = data.frame(x = x, y = y, id = id);
				return(xy);
			}
		}
		return(data.frame(x = x, y = y, id = id));
	});
	xyT = do.call(rbind, xyT);
	return(xyT);
}

