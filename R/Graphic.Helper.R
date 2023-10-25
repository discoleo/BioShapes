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

### Basic Math Functions:
# - reflect;
# - shift;


#' @export
as.bioshape = function(x) {
  if(inherits(x, "bioshape")) return(x);
  class(x) = c("bioshape", class(x));
  invisible(x);
}


################
### Geometry ###

### Distance
#' @export
dist.xy = function(x, y, as.sqrt = TRUE) {
  xy = (x[1] - x[2])^2 + (y[1] - y[2])^2;
  if(as.sqrt) xy = sqrt(xy);
  return(xy);
}


### Slope

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


### Basic Operations

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

### Translate All Points
# non-recursive list;
#' @export
shift.listNR = function(x, by = c(0, 0)) {
	for(i in seq_along(x)) {
			tmp = x[[i]];
			tmp$x = tmp$x + by[1];
			tmp$y = tmp$y + by[2];
			x[[i]] = tmp;
	}
	return(x);
}

### Shift Point or Line
# Shifts a point orthogonal to a given line;
# d = distance to shift (translate);
# simplify = TRUE: if(length(d) == 1) return a simple list(x, y);
#    otherwise return a data.frame;
# Note: direction of shift is not normalized with respect to quadrant;
#' @export
shiftLine = function(x, y, d = 1, slope = NULL,
		scale = 1, id.offset = 0, simplify = TRUE) {
	shift.ortho(x=x, y=y, d=d, slope=slope, scale=scale,
		id.offset=id.offset, simplify=simplify);
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
		scale=1, id.offset = 0, simplify = TRUE) {
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
    if(length(d) == 1 && simplify) {
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
    if(length(d) == 1 && simplify) {
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
  # Note: inconsistency when d = 1 between slope = Inf vs generic val;
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


######################

### Intersect Function

# Intersection between 2 line segments
# Returns:
#  - (x, y) of intersection(s);
#  - corresponding t1, t2, and
#  - n = number of (virtual) intersection points (0, 1, or 2 if overlapping);
# Real intersection only if: 0 <= t1 <= 1 & 0 <= t2 <= 1;
#' @export
intersect.lines = function(x, y, xB, yB) {
	xA1 = x[1]; xA2 = x[2]; xB1 = xB[1]; xB2 = xB[2];
	yA1 = y[1]; yA2 = y[2]; yB1 = yB[1]; yB2 = yB[2];
	#
	dyB = yB1 - yB2; dxB = xB1 - xB2;
	dyA = yA1 - yA2; dxA = xA1 - xA2;
	#
	div = dyB*dxA - dyA*dxB;
	# TODO: remaining special cases;
	if(div == 0) {
		# Parallel lines
		d = (xA1*yA2 - xA2*yA1)*dxB - (xB1*yB2 - xB2*yB1)*dxA;
		if(d == 0) {
			# Overlapping lines
			if(dxA == 0) {
				t1 = c(yB1 - yA2, yB2 - yA2) / dyA;
			} else {
				t1 = c(xB1 - xA2, xB2 - xA2) / dxA;
			}
			if(t1[1] > t1[2]) { tmp = t1[1]; t1[1] = t1[2]; t1[2] = tmp; }
			x = NA; y = NA; t1c = NA;
			if(t1[1] < 0) {
				# t1[1] = after start of 2nd segment;
				if(t1[2] >= 0) {
					t1c = c(0, min(1, t1[2]));
				}
			} else if(t1[1] <= 1) {
				t1c = c(t1[1], min(1, t1[2]));
			}
			if( ! is.na(t1c[1])) {
				x = xA1*t1c + (1-t1c)*xA2;
				y = yA1*t1c + (1-t1c)*yA2;
			}
			return(list(x=x, y=y, t1=t1, t2=Inf, n=2));
		} else {
			return(list(x=NA, y=NA, t1=Inf, t2=Inf, n=0));
		}
	}
	# Computations:
	t1  = (yA2 - yB1)*xB2 + dyB*xA2 + (yB2 - yA2)*xB1;
	t1  = - t1 / div;
	#
	div = dyB;
	t2  = yA2 - yB2 + dyA*t1;
	# TODO: div == 0
	t2  = t2 / div;
	#
	x = xA1*t1 + (1-t1)*xA2;
	y = yA1*t1 + (1-t1)*yA2;
	return(list(x=x, y=y, t1=t1, t2=t2, n=1));
}

#' @export
is.intersect.lines = function(x) {
	any(x$t1 >=0 & x$t1 <= 1 & x$t2 >= 0 & x$t2 <= 1);
}


###################

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

