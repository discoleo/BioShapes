###################
#
# Bachelor Thesis
#
# Title: BioShapes
#
# Candidate: Adrian Cotoc
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#
# in collaboration with Syonic SRL
# continous the work of Darian Voda
#
# GitHub: https://github.com/Adi131313/BioShapes

### Initial Ideas


####################

### Helper Functions

# - reflect;
# - shift;
# - plot;

#' @export
as.bioshape = function(x) {
  if(inherits(x, "bioshape")) return(x);
  class(x) = c("bioshape", class(x));
  invisible(x);
}

#' @export
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

### Shift Point or Line
# Shifts a point orthogonal to a given line;
# d = distance to shift (translate);
#' @export
shift.ortho = function(x, y, d = 1, slope = NULL,
			scale = 1, id.offset = 0) {
	shiftLine(x=x, y=y, d=d, slope=slope, scale=scale, id.offset=id.offset);
}
#' @export
shiftLine = function(x, y, d=1, slope=NULL,
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
      r = lapply(seq(along=d), function(id) {
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
      r = lapply(seq(along=d), function(id) {
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
  rez = lapply(seq(length(d)), function(id) shift.f(x, y, id))
  rez = data.frame(do.call(rbind, rez));
  return(rez);
}

# shift point p along line defined by (x, y);
# d = distance;
# TODO: Scale reference Ox
#' @export
shiftPoint = function(p, x, y, d=1, slope=NULL, scale=1) {
  if(is.null(slope)) {
    if(length(x) < 2 || length(y) < 2)
      stop("The base-line requires 2 points!");
    # TODO: handle if more than 2 points!
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
#' @export
is.helix.rev = function(phi, debug=FALSE) {
	# Warning: phi needs to be normalized to [0, 2*pi);
	phi.c = as.sin.intersect(phi[1]);
	isInc = phi.c[1] <= phi.c[2];
	if(debug) print(cbind(phi, phi.c));
	if(isInc) {
		# Note: strictly "<", but ">=";
		if(phi[2] < phi.c[1] || phi[2] >= phi.c[2]) {
			return(TRUE);
		}
	} else if(phi[2] > phi.c[2] && phi[2] <= phi.c[1]) {
		return(TRUE);
	}
	return(FALSE);
}


#######################

