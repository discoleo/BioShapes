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


### Functions to Generate Arrows

#####################

### Helper Functions

### Arrow function:
#' @export
arrow = function(x, y, type = "Simple", d=1, lwd=1, ...) {
  call = match.call();
  idType = match("type", names(call));
  if(is.na(idType)) {
    type = "Simple";
  } else {
    types = c("Simple", "Double", "Diamond",
              "Square", "Inverted", "T", "X", "N",
              "DoubleInverted", "Circle", "SolidSquare"); # TODO
    type = call[[idType]];
    type = pmatch(type, types);
    if(is.na(type)) stop("Invalid type!");
    call = call[ - idType];
    type = types[type];
  }
  # Function name:
  type = paste0("arrow", type);
  type = as.symbol(type);
  call[[1]] = type;
  eval(call)
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

#######################
### Specific Arrows ###

### Parameters:
# dH = horizontal length of ">";
# dV = vertical height of ">";
# d  = distance between each of ">>";
#      also equivalent to dH in the Simple types;
# d.lines = distance between multiple lines in the Tail;
# join = where to join Head & Tail;


### Simple Head
# - for consistency: join = 0;
#' @export
arrowSimple = function(x, y, d=-0.5, lwd=1, d.head=c(-d,d), d.lines=0,
                       h.lwd=lwd, col="red", scale=1, join=0, plot = TRUE) {
  slope = slope(x, y);
  qd = which.quadrant(x, y);
  sg = if(qd == 1 || qd == 4) 1 else -1;
  d  = sg * d;
  ### Head
  arrHead = list(H = arrowHeadSimple(x[2], y[2], slope=slope, d=d, dV = d.head, scale=scale),
               lwd = h.lwd);

  ### ArrowTail
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  if(any(d.lines != 0)) {
	arrow$xy = intersect.arrow(arrow$xy, arrHead[[1]]);
  }
  ### Full Arrow
  lst = list(Arrow=arrow, Head=arrHead, col=col);
  class(lst) = c("arrow", "list");
  # Plot lines:
  if(plot) lines(lst);
  invisible(lst);
}

### Double Lined Head: --->>
#' @export
arrowDouble = function(x, y, d=-0.5, lwd=1, d.head=-1, dV=c(-d.head, d.head), d.lines=0,
                       h.lwd=lwd, col="red", scale=1, join=0) {
  if(join > 2) stop("Unsupported value for join!");
  slope = slope(x, y);
  qd = which.quadrant(x, y);
  sg = if(qd == 1 || qd == 4) 1 else -1;
  d  = sg * d; d.head = sg * d.head;
  ### Head
  arrHead = arrowHeadDouble(x[2], y[2], slope=slope, d=d, dH=d.head, dV=dV, scale=scale);
  arrHead$lwd = h.lwd;
  ### ArrowTail
  if(join <= 1) {
    x[2] = arrHead[[2]]$x[2];
    y[2] = arrHead[[2]]$y[2];
  }
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope);
  if(any(d.lines != 0)) {
		# id vs ("H1", "H2"):
		id = if(join <= 1) 2 else 1;
		arrow$xy = intersect.arrow(arrow$xy, arrHead[[id]]);
  }
  ### Full Arrow
  lst = list(Arrow=arrow, Head=arrHead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}


### Inverted Head: ---<

### Simple Inverted Head
# if d <= 0: Head is within arrow-boundary;
#' @export
arrowInverted = function(x, y, d=-1, lwd=1, d.head=c(-d,d),
                         d.lines=0, h.lwd=lwd, col="red", scale=1, join=0) {
  slope = slope(x, y);
  isOK  = is.quadrant.right(x, y);
  ds = if(isOK) d else -d;
  tp = (d <= 0); # within boundary;
  ### Head
  arrHead = arrowHeadInverted(x[2], y[2], slope=slope, d=ds, dV=d.head,
		isTip = tp, scale=scale);
  arrHead = list(H = arrHead, lwd=h.lwd);
  ### Arrow Tail
  if(tp) {
    p = arrHead[[1]];
    x[2] = p$x[[2]];
    y[2] = p$y[[2]];
  }
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  if(any(d.lines != 0)) {
    # type = 3; # type = "Out";
	arrow$xy = intersect.arrow(arrow$xy, arrHead[[1]], type = 3);
  }
  ### Full Arrow
  lst = list(Arrow=arrow, Head=arrHead);
  class(lst) = c("arrow", "list");

  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}

### Double Lined Inverted Head
# dH = abs(d) ensures always inverted!
#' @export
arrowDoubleInverted = function(x, y, d=-0.25, lwd=1, dH=0.5, d.head=c(-dH, dH), d.lines=0,
                               h.lwd=lwd, col="red", scale=1, join=0) {
  if(join > 2) stop("Unsupported value for join!");
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  d  = if(isRQ) d else - d;
  ### Head
  arrHead = arrowHeadDoubleInverted(x[2], y[2], slope=slope, d=d, dH=dH, dV=d.head, scale=scale);
  midpoint = attr(arrHead, "Mid")
  arrHead$lwd = h.lwd;
  ### ArrowTail
  if(join <= 1) {
    midpoint = midpoint[[1]];
  } else {
    midpoint = midpoint[[2]];
  }
  x[2] = midpoint[1]
  y[2] = midpoint[2]
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  if(any(d.lines != 0)) {
		# type = 3; # type = "Out";
		id = if(join <= 1) 1 else 2;
		arrow$xy = intersect.arrow(arrow$xy, arrHead[[id]], type = 3);
  }
  ### Full Arrow
  lst = list(Arrow=arrow, Head=arrHead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}


### Multiple Lined Head: -->>>
# n = number of sub-components;
# d = distance between each "> >";
# dH = horizontal shift (shiftPoint)
# dV = vertical shift (shiftLine) of each ">";
#' @export
arrowN = function(x, y, n=1, d=-0.5, lwd=1, h.lwd=lwd, d.head=c(-d, d), d.lines=0,
                  col="red", scale=1, join=0) {
  if(join > n) stop("Unsupported value for join!");
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  d = if(isRQ) d else - d;
  ### Head
  arrHead = arrowHeadN(x[2], y[2], slope=slope, n=n, d=d, dV=d.head, scale=scale);
  arrHead$lwd = h.lwd;
  ### ArrowTail
  if(join == 0) {
    join = 1;
  }
  x[2] = arrHead[[join]]$x[2];
  y[2] = arrHead[[join]]$y[2];
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  if(any(d.lines != 0)) {
	arrow$xy = intersect.arrow(arrow$xy, arrHead[[join]]);
  }
  ### Full Arrow
  lst = list(Arrow=arrow, Head=arrHead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}


### Arrow w T-Head: ---|
#' @export
arrowT = function(x, y, d=0.2, lwd=1, d.head=c(-d, d), d.lines=0, h.lwd=lwd,
                  col="red", scale=1, join=0, lty=1) {
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  d = if(isRQ) d else - d;
  ### Head
  if(is.list(d.head)) {
    ahead = lapply(seq_along(d.head), function(id) {
      arrowHeadT(x[2], y[2], slope=slope, dV=d.head[[id]], scale=scale)
    });
    ahead$lwd = h.lwd;
  } else {
    ahead = list(arrowHeadT(x[2], y[2], slope=slope, dV=d.head, scale=scale),
                 lwd = h.lwd);
  }
  ### ArrowTail
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  ### Full Arrow
  lst = list(Arrow=arrow, Head=ahead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col, lty=lty);
  invisible(lst);
}

### Measurements

### Arrow for Measurements
#' @export
arrowMeasure = function(x, y, d=-0.5, lwd=1, d.head=c(-d,d), dT=d.head, d.lines=0,
                        h.lwd=lwd, col="red", scale=1, join=0) {
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  d = if(isRQ) d else - d;
  ### Head
  arrHead = arrowHeadMeasure(x[2], y[2], slope=slope, d=d, dV = d.head, dT=dT, scale=scale);
  arrHead$lwd = h.lwd;
  ### ArrowTail
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  if(any(d.lines != 0)) {
	arrow$xy = intersect.arrow(arrow$xy, arrHead[[1]]);
  }
  ### Full Arrow
  lst = list(Arrow=arrow, Head=arrHead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}

### Measure: |<--->|
# - Double Headed Arrow;
#' @export
measure = function(x, y, type=c("in", "out"), lty=1, lwd=1, col=1,
                   d=c(-1,1), dH=0.5, d.head=c(-dH, dH), d.lines=0,
                   lwd.head=lwd, scale=1, plot = TRUE) {
  type  = match.arg(type);
  slope = slope(x, y);
  # Head
  if(type == "in") {
    h1 = arrowHeadSimple(x[1], y[1], slope=slope, d = abs(dH), dV = - d.head, scale=scale);
    h2 = arrowHeadSimple(x[2], y[2], slope=slope, d = - abs(dH), dV = d.head, scale=scale);
  } else {
    h1 = arrowHeadSimple(x[1], y[1], slope=slope, d = - abs(dH), dV = - d.head, scale=scale);
    h2 = arrowHeadInverted(x[2], y[2], slope=slope, d = abs(dH), dV = d.head, scale=scale);
  }
  t1 = arrowHeadT(x[1], y[1], slope=slope, dV = d, scale=scale);
  t2 = arrowHeadT(x[2], y[2], slope=slope, dV = d, scale=scale);
  hh = list(h1, h2, t1, t2, lwd=lwd.head);
  # Tail
  arrTail = arrowTail(x, y, slope=slope, d.lines=d.lines, lwd=lwd);
  # arrTail$lty = lty; # TODO
  # Full Arrow
  lst = list(Tail=arrTail, Head=hh);
  class(lst) = c("arrow", "list");
  if(plot) {
    lines.arrow(lst, col=col, lty=lty);
  }
  invisible(lst)
}


### Other Head Types

### Arrow X: ---X
# - for consistency: join = 0;
#' @export
arrowX = function(x, y, d=0.5, lwd=1, d.head=c(-d, d), d.lines=0,
                  h.lwd=lwd, col="red", scale=1, join=0) {
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  d = if(isRQ) d else - d;
  ### Head
  arrHead = arrowHeadX(x[2], y[2], slope=slope, d = - d, dV = d.head, scale=scale);
  ahead   = list(H = arrHead, lwd = h.lwd);
  midHxy  = attr(arrHead, "Mid");
  ### ArrowTail
  x[2] = midHxy[1];
  y[2] = midHxy[2];
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  if(any(d.lines != 0)) {
		isSimple = nrow(arrHead) == 4;
		xyH = if(isSimple) arrHead
			else list(x = arrHead$x[c(4,5,6)], y = arrHead$y[c(4,5,6)]);
		arrow$xy = intersect.arrow(arrow$xy, xyH);
  }
  ### Full Arrow
  lst = list(Arrow=arrow, Head=ahead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}


### Arrow Circle: ---O
#' @export
arrowCircle = function(x, y, r=0.5, lwd=1, d.lines=0,
                       h.lwd=lwd, col="red", fill=NULL, scale=1, join=0) {
  if(join > 3) stop("Unsupported value for join!");
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  r = if(isRQ) r else - r;
  ### Head
  arrHead = arrowHeadCircle(x[2], y[2], slope=slope, r=r, scale=scale);
  mid = attr(arrHead, "start")
  arrHead$lwd = h.lwd;
  arrHead[[1]]$fill = fill;
  ### ArrowTail
  if(join < 3) {
    mid = if(join == 2) mid[[2]] else mid[[1]];
    x[2] = mid[1];
    y[2] = mid[2];
  }
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  ### Full Arrow
  lst = list(Arrow=arrow, Head=arrHead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}

### Arrow Diamond: ---<>
#' @export
arrowDiamond = function(x, y, d=0.2, lwd=1, d.head=c(-1, 1), d.lines=0,
		h.lwd=lwd, col="red", fill=NULL, scale=1, join=0) {
  if(join > 2) stop("Unsupported value for join!");
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  d = if(isRQ) d else - d;
  ### Head
  arrHead = list(
		H = arrowHeadDiamond(x[2], y[2], slope=slope, d=d, dV=d.head, scale=scale),
		lwd = h.lwd);
  ### ArrowTail
  if(join == 0 || join == 1) {
    x[2] = arrHead[[1]]$x[2];
    y[2] = arrHead[[1]]$y[2];
  }
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  if(any(d.lines != 0)) {
    # type = 3; # type = "Out";
	arrow$xy = intersect.arrow(arrow$xy, arrHead[[1]], type = 3);
  }
  ### Full Arrow
  if( ! is.null(fill)) {
		xy   = arrHead$H;
		xy$x = xy$x[ - 5];
		xy$y = xy$y[ - 5];
		xy$fill = fill;
		class(xy) = c("polygon", "list");
		arrHead$H = xy;
  }
  lst = list(Arrow = arrow, Head = arrHead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}


### Arrow Square: ---[]
# default: fill = NULL;
#' @export
arrowSquare = function(x, y, d=-0.5, lwd=1, d.head=c(d, -d)/2, d.lines=0,
                       h.lwd=lwd, col="red", fill=NULL, scale=1, join=0) {
  if(join > 2) stop("Unsupported value for join!");
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  isLow = (d < 0);
  d = if(isRQ) d else - d;
  ### Head
  arrHead = arrowHeadSquare(x[2], y[2], slope=slope, d=d, dV=d.head, scale=scale);
  if( ! is.null(fill)) {
    arrHead$fill = fill;
    class(arrHead) = c("polygon", class(arrHead));
  }
  ahead = list(H = arrHead, lwd = h.lwd);
  ### ArrowTail
  if(join < 2 && isLow) {
    mid  = attr(arrHead, "Mid");
    x[2] = mid[1]; y[2] = mid[2];
  }
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  ### Full Arrow
  lst = list(Arrow=arrow, Head=ahead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}

### Arrow Solid Square
#' @export
arrowSolidSquare = function(x, y, d=-0.5, lwd=1, d.head=c(d, -d)/2, d.lines=0,
                            h.lwd=lwd, col="red", fill=col, scale=1, join=0) {
  if(join > 2) stop("Unsupported value for join!");
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  isLow = (d < 0);
  d = if(isRQ) d else - d;
  ### Head
  arrHead = arrowHeadSquare(x[2], y[2], slope=slope, d=d, dV=d.head, scale=scale);
  arrHead$fill = fill;
  class(arrHead) = c("polygon", class(arrHead));
  ahead = list(H = arrHead, lwd = h.lwd);
  ### ArrowTail
  if(join < 2) {
    mid  = attr(arrHead, "Mid");
    x[2] = mid[1]; y[2] = mid[2];
  }
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  ### Full Arrow
  lst = list(Arrow=arrow, Head=ahead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}


### Arrow Triangle: ---|>
#' @export
arrowTriangle = function(x, y, d=-0.5, lwd=1, d.head=c(-d,d), d.lines=0,
                         h.lwd=lwd, col="red", fill = NULL, scale=1, join=0) {
  if(join > 2) stop("Unsupported value for join!");
  slope = slope(x, y);
  isRQ  = is.quadrant.right(x, y);
  isLow = (d < 0);
  d = if(isRQ) d else - d;
  ### Head
  arrHead = arrowHeadTriangle(x[2], y[2], slope=slope, d=d, dV = d.head, scale=scale);
  mid     = attr(arrHead, "Mid");
  arrHead = list(H = arrHead, lwd = h.lwd);
  if( ! is.null(fill)) {
	arrHead[[1]]$fill = fill;
	class(arrHead[[1]]) = c("polygon", "list");
  }
  ### ArrowTail
  if(join < 2 && isLow) {
    x[2] = mid[1]; y[2] = mid[2];
  }
  arrow = arrowTail(x, y, d.lines=d.lines, lwd=lwd, slope=slope, scale=scale);
  ### Full Arrow
  lst = list(Arrow=arrow, Head=arrHead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}

#################

### Arrow: Tail = Square Wave
#' @export
arrowSquareWave = function(x, y, n, d=-0.5, dV=c(-1,1), d.head=c(-d,d),
                           col=1, lwd=1, h.lwd=lwd, scale=1, join=0) {
  if(join > 2) stop("Unsupported value for join!");
  slope = slope(x, y);
  ### Head
  ahead = list(arrowHeadSimple(x[2], y[2], slope=slope, d=d, dV=d.head, scale=scale),
               lwd = h.lwd);
  ### ArrowTail
  l = shiftLine(x, y, d=dV, slope=slope, scale=scale);
  # Split lines:
  p1 = c(x[1], y[1]);
  p2 = c(x[2], y[2]);
  # l0 = split.line(p1, p2, n+1);
  l1 = l[l$id == 1, c("x", "y")];
  l1 = split.line(l1[,1], l1[,2], n+1);
  l2 = l[l$id == 2, c("x", "y")];
  l2 = split.line(l2[,1], l2[,2], n+1);
  # Join lines:
  join.lines = function(l, l1, l2) {
    pos = length(l1) - 1; # before last
    mid = (l1[pos] + l2[pos]) / 2;
    if(n %% 2 == 1) {
      l1 = head(l1, -1);
      l2 = c(tail(head(l2, -2), -1), mid, mid);
    } else {
      l1 = head(l1, -2);
      l2 = c(tail(head(l2, -1), -1));
    }
    l1 = matrix(l1, nrow=2);
    l2 = matrix(l2, nrow=2);
    lt = rbind(l1, l2);
    lt = c(as.vector(lt));
    if(n %% 2 == 0) lt = c(lt, mid);
    lt = c(l[1], lt, l[2]);
  }
  x = join.lines(x, l1$x, l2$x);
  y = join.lines(y, l1$y, l2$y);
  arrow = list(list(x=x, y=y), lwd=lwd);
  ### Full Arrow
  lst = list(Arrow=arrow, Head=ahead);
  class(lst) = c("arrow", "list");
  # Plot lines:
  lines(lst, col=col);
  invisible(lst);
}


######################

### Specialized Arrows

### Chemical Reactions: Double Halves
#' @export
arrowDHalf = function(x, y, d = 0.1, dH = 0.5, d.head = c(0, dH), scale = 1,
		lwd = 1, col = NULL, plot = FALSE) {
	slope = slope(x, y);
	if(length(d) == 1) d = c(d, -d);
	qd = which.quadrant(x, y);
	l1 = shiftLine(x, y, d = d[1], slope=slope, scale=scale);
	l2 = shiftLine(rev(x), rev(y), d = d[2], slope=slope, scale=scale);
	#
	if(length(lwd) == 1) lwd = c(lwd, lwd);
	if(length(col) == 1) col = c(col, col);
	a1 = arrowSimple(l1$x, l1$y, d = -dH, d.head = d.head, col=col[[1]], lwd=lwd[[1]],
		scale=scale, plot = FALSE);
	a2 = arrowSimple(l2$x, l2$y, d = -dH, d.head = - d.head, col=col[[2]], lwd=lwd[[2]],
		scale=scale, plot = FALSE);
	lst = list(A1 = a1, A2 = a2);
	class(lst) = c("multiArrow", "list");
	# Plot arrow:
	if(plot) lines(lst);
	invisible(lst);
}




### Generate complex lines

### Complex Lines:
# - just a simple Example;
# TODO: correct function
lineBanded = function(x, y, w=0.1, delta=0.25, lwd=1.5, lty=1, n=NULL, col="black", slope=NULL) {
  if(is.null(slope)) {
    slope = slope(x, y);
  }
  lsh = shiftLine(x, y, d=w, slope=slope);
  distxy = sqrt((x[1] - x[2])^2 + (y[1] - y[2])^2);
  if(is.null(n)) {
    # TODO: use also lwd for better accuracy;
    n = distxy / delta;
  }
  # partition lines:
  dx = lsh[3,1] - lsh[1,1]; dx = dx / n;
  dy = lsh[3,2] - lsh[1,2]; dy = dy / n;
  xup = c(lsh[1,1], lsh[1,1] + dx*seq(n-1), lsh[3,1]);
  xdn = c(lsh[2,1], lsh[2,1] + dx*seq(n-1), lsh[4,1]);
  yup = c(lsh[1,2], lsh[1,2] + dy*seq(n-1), lsh[3,2]);
  ydn = c(lsh[2,2], lsh[2,2] + dy*seq(n-1), lsh[4,2]);
  # Plot:
  for(id in seq(n)) {
    lines(c(xup[id], xdn[id]), c(yup[id], ydn[id]), lwd=lwd, lty=lty, col=col);
  }
  invisible(cbind(xup, xdn, yup, ydn));
}

### Braces ###
#' @export
braces.curly = function(center = c(0,0), scale = c(6, 1), left.open = TRUE, th = 90,
                        pow=2.5, mid = 0.625, y0.scale = 1.75, npx = 32) {
  x0 = seq(-1, 0, length.out = npx);
  dx = x0 + mid;
  sg = sign(dx);
  if( ! left.open) sg = - sg;
  y = sg * abs(dx)^pow;
  y[1] = y0.scale * y[1];
  y = c(y, rev(y)) * scale[2];
  len = scale[1] / 2;
  x = seq(- len, len, length.out = 2*npx);
  if(th != 0) {
    th = th * pi / 180;
    sn = sin(th); cs = cos(th);
    r = cbind(x*cos(th) - y*sin(th), x*sin(th) + y*cos(th));
  } else {
    r = cbind(x, y);
  }
  r[,1] = r[,1] + center[1];
  r[,2] = r[,2] + center[2];
  return(r);
}

### Braces 2 ###
#' @export
braces.curly2 = function(center = c(0, 0), cex=4, left.open = TRUE, th=90, ...) {
  txt = if(left.open) "{" else "}";
  text(x = center[1], y = center[2], txt, cex=cex, srt = th - 90, ...);
}
