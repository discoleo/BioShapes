#######################################
#
# Title: BioShapes
#
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
# GitHub: https://github.com/Adi131313/BioShapes
#
# 2. Bachelor Thesis: Darian Voda (2021-2022)



####################

### Helper Functions

### Basic Math Functions:
# - reflect;
# - shift;


#' @export
as.bioshape = function(x) {
  if(inherits(x, "bioshape")) return(x);
  class(x) = c("bioshape", class(x));
  invisible(x);
}


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


### Negate Formula
#' @export
formula.neg = function(x) {
	tmp = expression(- (0))[[1]];
	tmp[[2]][[2]] = x[[2]];
	x[[2]] = tmp;
	return(x)
}


################
### Geometry ###

# Deprecated:
compute_slope = function(x, y) {
  warning("Deprecated");
  return(slope(x, y));
}

#' @export
slope = function(x, y) {
  if(length(x) < 2 || length(y) < 2) {
    stop("The base-line requires 2 points!");
    # TODO: check if it works properly.
  }

  if(length(x) == 2){
    slope = (y[[2]] - y[[1]]) / (x[[2]] - x[[1]]);
  }
  else{
    slope = diff(y) / diff(x);
  }
  return(slope);
}

### Quadrant
# 2|1
# 3|4
#' @export
which.quadrant = function(x, y) {
	qdr = if(x[1] <= x[2]) {
		if(y[1] <= y[2]) 1 else 4;
	} else {
		if(y[1] < y[2]) 2 else 3;
	}
	return(qdr);
}
#' @export
which.quadrant.phi = function(phi) {
	if(phi < 0 || phi > 2*pi) phi = as.radians0(phi);
	if(phi <= pi/2) return(1);
	if(phi < pi) return(2);
	if(phi <= 3*pi/2) return(3);
	return(4);
}


### Center: 4 points
# t = partition;
#' @export
center.p4 = function(p1, p2, p3, p4, t = c(1/2, 1/2)) {
	tx = cbind(t[1], 1 - t[1]);
	ty = rbind(t[2], 1 - t[2]);
	xx = matrix(c(p1[1], p2[1], p3[1], p4[1]), 2);
	yy = matrix(c(p1[2], p2[2], p3[2], p4[2]), 2);
	mid.x = as.vector((tx %*% xx) %*% ty);
	mid.y = as.vector((tx %*% yy) %*% ty);
	center = c(mid.x, mid.y);
	return(center);
}


### Reflect Point across line
# p = c(x, y) # the point to reflect;
#' @export
reflect = function(x, y, p, slope=NULL) {
  if(is.null(slope)) {
    if(length(x) < 2 || length(y) < 2)
      stop("The base-line requires 2 points!");
    # TODO: handle if more than 2 points!
    slope = slope(x,y) # (y[[2]] - y[[1]]) / (x[[2]] - x[[1]]);
  }

  if(slope == 0) {
    # Horizontal Line
    return(c(p[1], 2*y[1]-p[2]))
  } else if(x[[1]] == x[[2]]) {
    # Vertical Line
    return(c(2*x[1]-p[1], p[2])) # De ce?
    # TODO: distanta de la punct la dreapta + punctul initial ?
  }

  sl.orto = - 1 / slope;
  # intersection point:
  div = slope - sl.orto; # div always > 0!
  x.int = (x[1]*slope - p[1]*sl.orto - y[1] + p[2]) / div;
  y.int = (x.int - p[1])*sl.orto + p[2];
  print(c(x.int, y.int))
  # Reflected point:
  x.rfl = 2*x.int - p[1];
  y.rfl = 2*y.int - p[2];
  return(c(x.rfl, y.rfl));
}


### Shift Points

#' @export
shift = function(x, y, d = 1, slope = NULL, scale = 1, ...) {
	UseMethod("shift");
}
#' @export
shift.point = function(p, x, y, d = 1, slope = NULL, scale = 1, ...) {
	UseMethod("shift.point");
}

### Shift Point or Line
# Shifts a point orthogonal to a given line;
# d = distance to shift (translate);
# Note: direction of shift is not normalized with respect to quadrant;
#' @export
shiftLine = function(x, y, d = 1, slope = NULL,
		scale = 1, id.offset = 0) {
	shift.ortho(x=x, y=y, d=d, slope=slope, scale=scale, id.offset=id.offset);
}
#' @export
shift.ortho.df = function(xy, d = 1, slope = NULL,
		scale = 1, id.offset = 0) {
	xy = shift.ortho(xy$x, xy$y, d=d, slope=slope,
		scale=scale, id.offset=id.offset);
	return(xy);
}
#' @export
shift.ortho = function(x, y, d=1, slope=NULL,
		scale=1, id.offset = 0) {
  if(is.null(slope)) {
    if(length(x) < 2 || length(y) < 2)
      stop("The base-line requires 2 points!");
    # TODO: handle if more than 2 points!
    slope = slope(x,y);
  } else {
    if(missing(y)) {
      # both coordinates encoded using parameter x;
      y = x[2]; x = x[1];
    }
  }
  ### Vertical Line
  if(abs(slope) == Inf) {
    if(length(d) == 1) {
      r = data.frame(x = x - d, y = y);
    } else {
      r = lapply(seq(along = d), function(id) {
        data.frame(x = x - d[id], y = y, id = id + id.offset);
      })
      r = do.call(rbind, r);
    }
    return(r)
  }
  ### Horizontal Line
  if(slope == 0) {
    d = d * scale;
    if(length(d) == 1) {
      r = data.frame(x = x, y = y + d);
    } else {
      r = lapply(seq(along = d), function(id) {
        data.frame(x = x, y = y + d[id], id = id + id.offset);
      })
      r = do.call(rbind, r);
    }
    return(r)
  }
  ### Oblique Line
  sl.orto = - 1 / slope;
  sl2 = sqrt(sl.orto^2 + 1/scale^2);
  sl.orto = sl.orto * scale;
  # shift Start- & End-points:
  shift.f = function(x, y, id) {
    delta = d[id] / sl2;
    x.sh  = x - delta / scale;
    y.sh  = y - delta*sl.orto;
    data.frame(x=x.sh, y=y.sh, id = id + id.offset);
  }
  rez = lapply(seq(along = d), function(id) shift.f(x, y, id))
  rez = data.frame(do.call(rbind, rez));
  return(rez);
}

# shift point p along line defined by (x, y);
# d = distance;
# TODO: Scale reference Ox
#' @export
shift.points.df = function(p, x, y, d=1, slope=NULL, scale=1, simplify = TRUE) {	
	len = nrow(p);
	if(len == 0) return(list());
	if(is.null(slope)) {
		if(length(x) < 2 || length(y) < 2)
			stop("The base-line requires 2 points!");
		slope = slope(x,y);
	}
	lst = lapply(seq(len), function(id) {
		pp = shiftPoint(p[id, c("x", "y")], slope=slope, d=d, scale=scale);
		data.frame(pp);
	})
	if( ! simplify) return(lst);
	lst = do.call(rbind, lst);
	names(lst) = c("x", "y");
	if(match("id", names(p))) {
		lst$id = rep(p$id, each = length(d));
	}
	return(lst);
}
shiftPoint = function(p, x, y, d = 1, slope = NULL, scale = 1) {
	shift.point(p=p, x=x, y=y, d=d, slope=slope, scale=scale);
}
#' @export
shift.point.default = function(p, x, y, d = 1, slope = NULL, scale = 1) {
  if(is.null(slope)) {
    if(length(x) < 2 || length(y) < 2)
      stop("The base-line requires 2 points!");
    slope = slope(x,y);
  }
  if(length(p) < 2) stop("Point needs both x & y coordinates!");
  # V line: if(x[1] == x[2])
  sgn = sign(slope);
  if(abs(slope) == Inf) {
    r = cbind(x = p[1], y = p[2] + d*scale*sgn);
    return(r)
  }
  slope.sqrt = scale / sqrt(slope^2 + scale^2);
  dx = d * slope.sqrt;
  dy = dx * slope;
  xsh = p[[1]] + dx;
  ysh = p[[2]] + dy;
  return(cbind(x=xsh, y=ysh));
}

#### Helper Functions for Circles ####
#' @export
solveCircleIntersection = function(mid1, mid2, r1, r2, digits=4, debug=TRUE) {
  xp = mid1[[1]]; yp = mid1[[2]];
  xc = mid2[[1]]; yc = mid2[[2]];
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
  if(debug) print(c(a, b1, b0, D));
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


### Split Line
#' @export
split.line = function(x, y, n) {
  if(n == 1) {
    lst = list(x=x, y=y);
    return(lst);
  }
  if(n <= 0) stop("Wrong number of fragments!");
  # Split
  t  = seq(n - 1) / n;
  rt = rev(t);
  xs = x[1]*rt + x[2]*t;
  ys = y[1]*rt + y[2]*t;
  xs = c(x[1], xs, x[2]);
  ys = c(y[1], ys, y[2]);
  lst = list(x=xs, y=ys);
  return(lst);
}

### Split into alternate & separate lines
#' @export
split.AltLines.matrix = function(xy) {
	len = nrow(xy);
	xyS = xy[seq(1, len, by=2), , drop=FALSE];
	xyE = xy[seq(2, len, by=2), , drop=FALSE];
	xyS = as.data.frame(xyS);
	xyE = as.data.frame(xyE);
	len = nrow(xyS);
	xyS$id = seq(len);
	xyE$id = seq(len);
	xy = rbind(xyS, xyE);
	names(xy) = c("x", "y", "id");
	return(xy);
}


### Distance
#' @export
dist.xy = function(x, y, as.sqrt = TRUE) {
  xy = (x[1] - x[2])^2 + (y[1] - y[2])^2;
  if(as.sqrt) xy = sqrt(xy);
  return(xy);
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

# p1 = Translate to p1;
#' @export
rotate = function(x, y, slope, p1=c(0,0)) {
  if(abs(slope) == Inf) {
    sgn = sign(slope);
    dx = y; dy = x * sgn;
    lst = list(x = p1[1] + dx, y = p1[2] + dy);
    return(lst);
  }
  sdiv = 1 / sqrt(slope^2 + 1);
  # Rotation matrix: by column
  # rotm = matrix(sdiv * c(1, s, -s, 1), ncol=2, nrow=2);
  dx = (x - slope*y) * sdiv; # + p1[1];
  dy = (slope*x + y) * sdiv; # + p1[2];
  lst = list(x = p1[1] + dx, y = p1[2] + dy);
}


### Helix / Spirals

### Radians: between [0, 2*pi);
#' @export
as.radians0 = function(x, tol=1E-10) {
	p2 = 2*pi;
	xf = floor(x / p2);
	x  = x - p2*xf;
	x[abs(x) < tol] = 0;
	return(x);
}

### Sinusoid Origins
# - which sinusoids intersect given sinusoid (same origin);
#' @export
as.sin.intersect = function(phi) {
	sapply(phi, function(phi) {
		as.radians0(pi - phi + c(0, pi));
	})
}

# Intersection of 2 shifted-Sine Functions
# from, to = integers (cycles of 2*pi);
#' @export
which.intersect.sin = function(phi, n, from = 0, to = NULL) {
  phi.eq = (pi*c(1,3) - sum(phi))/2;
  id.neg = which(phi.eq < 0);
  if(length(id.neg) > 0) {
    shift = abs(floor(phi.eq[id.neg] / (2*pi)));
    phi.eq[id.neg] = phi.eq[id.neg] + (2*pi)*shift;
  }
  # Sort:
  if(phi.eq[2] < phi.eq[1]) phi.eq = rev(phi.eq);
  # Last:
  phi.rad = phi.eq / (2*pi);
  if(is.null(to)) {
    # complete 2*pi cycles: starts at 0;
	# tol = 0.01;
    to = floor(n - phi.rad[2] + 0.01);
  }
  x0 = sapply(seq(from, to, by=1), function(n) {
    phi.eq + 2*n*pi;
  })
  last = n - 1 - to - phi.rad[1] + 0.01;
  # print(c(to, last, phi.rad));
  if(last >= 0)
    x0 = c(x0, x0[length(x0)] + pi);

  x0 = sort(x0);
  x1 = x0 + phi[1];
  # x0 = non-shifted; x1 = shifted;
  return(list(x0 = x0, x1=x1));
}

### Order of 2 Helices
# phi = phase shift of the 2 helices;
#' @export
is.helix.rev = function(phi, debug=FALSE) {
	# Warning: phi needs to be normalized to [0, 2*pi);
	# Note: Case phi[2] == phi[1] is indeterminate;
	if(phi[1] < pi/2) {
		isRev = (phi[2] > phi[1]) && (phi[2] <= pi - phi[1]);
		return(isRev);
	}
	if(phi[1] == pi/2) return(FALSE); # always FALSE;
	if(phi[1] < pi) {
		isRev = (phi[2] < phi[1]) && (phi[2] > pi - phi[1]);
		return(isRev);
	}
	if(phi[1] < 3*pi/2) {
		isRev = (phi[2] < phi[1]) || (phi[2] > 3*pi - phi[1]);
		return(isRev);
	}
	if(phi[1] == 3*pi/2) return(TRUE); # always TRUE;
	# Case: phi[1] > 3*pi/2;
	isRev = phi[2] > phi[1] || phi[2] <= 3*pi - phi[1];
	return(isRev);
	
	# TODO: [old] cleanup?
	# phi.c = as.sin.intersect(phi[1]);
	# isInc = phi.c[1] <= phi.c[2];
}


#######################

