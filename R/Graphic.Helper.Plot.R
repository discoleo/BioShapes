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

### Create new plot window:
# Convenience function:
#' @export
plot.base = function(xlim=c(-2,10), ylim=c(-2,10), axt=c(1,2), asp=1, set.par=TRUE,
		cex.axis = 1, ...) {
  hasAxis = ! is.null(axt);
  if(set.par) {
    mar = 0.25 + if(hasAxis) c(2,2,0,0) else c(0,0,0,0);
    par.old = par(mar = mar);
  } else par.old = NA;
  plot.new()
  plot.window(xlim=xlim, ylim=ylim, asp=asp, ...);
  if(hasAxis) {
    lapply(axt, function(a) axis(a, cex.axis = cex.axis));
  }
  invisible(par.old);
}

### Helper

#' @export
plot.circle = function(r, center=c(0,0), col=1, fill=NULL, N=128, ...) {
  x = r * cos(2*pi*seq(0, N-1) / N) + center[1];
  y = r * sin(2*pi*seq(0, N-1) / N) + center[2];
  polygon(x, y, border=col, col=fill, ...);
  invisible();
}

# Note: mainly for testing the hashed lines;
#' @export
plot.circle.hash = function(n, phi = pi/3, r = 1, center = c(4,4),
		col.line = "green", scale = 1, ...) {
	tmp = circle.hash(n, phi=phi, r=r, center=center, scale=scale, ...)
	plot.base(asp = 1/scale);
	lines(tmp);
	if( ! is.null(col.line)) {
		abline(v = center[1] + c(-r,r), col = col.line);
	}
}

#' @export
plot.circle.arc = function(r, center, phi, col=1, fill=NULL, ...) {
  shape::plotcircle(r, mid=center, from=phi[1], to=phi[2], lcol=col, col=fill, ...);
}

#' @export
plot.ellipse = function(r, center, phi = c(0, pi), theta = 0, lwd = 1,
		col = 1, fill = NULL, scale = 1, ..., N = 64) {
	# shape::plotellipse(rx = r[1], ry = r[2], mid = center, angle = th * 180 / pi,
	#	from = phi[1], to = phi[2], lwd=lwd, lcol = col, col = fill, ...);
	dr = max(r) / N;
	xy = shape::getellipse(rx = r[1], ry = r[2], mid = center, angle = theta * 180 / pi,
		from = phi[1], to = phi[2], dr = dr);
	if(scale != 1) xy = shape::rotatexy(xy, angle = 0, mid = center, asp = TRUE);
	if ( ! is.null(fill))
        polygon(xy, col = fill, border = NA);
	lines(xy, lwd=lwd, col=col, ...);
}

### Plot:

### Bio-Shapes
#' @export
lines.bioshape = function(x, lwd=NULL, col=NULL, ...) {
  lines.object.base(x, lwd=lwd, col=col, ...)
  invisible();
}

# Note: represents usually an error;
#' @export
lines.list = function(x, y, lwd=NULL, ...) {
  warning("Class is missing: Only list !");
  sapply(x, function(lst) {
    lwd = if(is.null(lwd)) {
      if(is.null(lst$lwd)) 1 else lst$lwd;
    } else lwd;
    lines(lst$x, lst$y, lwd=lwd, ...);
  });
}

### Base function
#' @export
lines.object.base = function(x, lwd, col, fill=NULL, ...) {
	# do NOT overwrite user-value;
	if(is.null(lwd)) {
		lwd = if(is.null(x$lwd)) 1 else x$lwd;
	}
	if(is.null(col)) {
		col = if(is.null(x$col)) 1 else x$col;
	}
	x$lwd = NULL; x$col = NULL;
	# Actual components:
	basef = function(lst, lwd, col, ...) {
		if(! is.null(lst$lwd)) { lwd = lst$lwd; lst$lwd = NULL; }
		if(! is.null(lst$col)) { col = lst$col; lst$col = NULL; }
		#
		if(inherits(lst, "bioshape")) {
			lapply(lst, basef, lwd = lwd, col = col, ...);
			return();
		}
		if(inherits(lst, "circle")) {
			if(is.null(fill)) fill = lst$fill;
			if(inherits(lst$center, "matrix")) {
				lapply(seq(nrow(lst$center)), function(nr) {
					shape::plotellipse(rx = lst$r, ry = lst$r, mid = lst$center[nr, ],
						lcol=col, col=fill, lwd=lwd, ...);
				})
			} else {
				shape::plotellipse(rx = lst$r, ry = lst$r,
					mid = lst$center, lcol=col, col=fill, lwd=lwd, ...);
			}
	} else if(inherits(lst, "ellipse")) {
		fill = lst$fill;
		if(inherits(lst$center, "matrix")) {
			print("TODO")
			lapply(seq(nrow(lst$center)), function(nr) {
				# TODO
			})
		} else {
			lst$col = col; lst$fill = fill;
			lst$lwd = lwd; # TODO: check effects;
			do.call("plot.ellipse", lst);
			# plot.ellipse(r = lst$r, center = lst$center, phi = lst$phi, th = lst$th,
			#	lwd=lwd, col=col, fill=fill, scale = lst$scale, ...);
		}
    } else if(inherits(lst, "circle.arc")) {
      if(is.null(fill)) fill = lst$fill;
      if(inherits(lst$center, "matrix")){
        lapply(seq(nrow(lst$center)), function(nr){
          plot.circle.arc(
            r = lst$r, center = lst$center[nr, ], phi = lst$phi[nr, ],
            col=col, fill=fill, lwd=lwd, ...);
        })
      } else {
        plot.circle.arc(
          r = lst$r, center = lst$center, phi = lst$phi,
          col=col, fill=fill, lwd=lwd, ...);
      }

    } else if(inherits(lst, "polygon")) {
      # if(is.null(fill)) fill = lst$fill;
      fill = if(is.null(lst$fill)) NA else lst$fill;
      polygon(lst$x, lst$y, col=fill, border=col, lwd = lwd, ...);
    } else if(inherits(lst, "text")) {
		len = length(lst);
		if(len == 3) { text(lst$x, lst$y, lst$labels, col=col); }
		else {
			if( ! is.null(col)) lst$col = col;
			do.call("text", lst);
		}
	} else if(inherits(lst, "arrow")) {
		lines.arrow(lst, col=col);
	} else {
		# warning("Only lines");
		if(is.null(lwd)) lwd = 1;
		# Note: lines.list uses the same color for all elements!
		if(inherits(lst, "lines.list")) {
			lapply(lst, function(lst) {
				tmp = lst; tmp$lwd = lwd; tmp$col = col;
				tmp = merge.first.list(lst, list(...));
				do.call("lines", tmp);
				# lines(lst$x, lst$y, lwd=lwd, col=col, ...);
			});
		} else if(inherits(lst, "data.frame")) {
			lst = split(lst[, c("x", "y")], lst$id);
			lapply(lst, function(lst) {
				lines(lst$x, lst$y, lwd=lwd, col=col, ...);
			});
		} else {
			lst$lwd = lwd; lst$col = col;
			lst = merge.first.list(lst, list(...));
			if(length(lst$col) == 1) {
				do.call("lines", lst);
			} else {
				len = length(lst$x);
				lapply(seq(len - 1), function(id) {
					id = c(id, id + 1);
					lst$x = lst$x[id]; lst$y = lst$y[id];
					lst$col = lst$col[id];
					do.call("lines", lst);
				});
			}
			# lines(lst$x, lst$y, lwd=lwd, col=col, ...);
		}
    }
  }
  lapply(x, basef, lwd=lwd, col=col, ...);
  invisible();
}

### list(Tail=list(...), Head=list(...))
#' @export
lines.arrow = function(x, lwd = NULL, col = NULL, lty = NULL, ...) {
	if(is.null(col)) {
		if(is.null(col <- x$col)) col = 1;
	}
	if(is.null(lty)) {
		if(is.null(lty <- x$lty)) lty = 1;
	}
	### ArrowTail
	arrow = x$Arrow;
	lines.object.base(arrow, lwd=lwd, col=col, lty=lty, ...)
	### ArrowHead
	ahead = x$Head;
	lines.object.base(ahead, lwd=lwd, col=col, ...);
	#
	invisible();
}
#' @export
lines.multiArrow = function(x, lwd=NULL, col = NULL, lty = 1, ...) {
	lapply(x, function(x) lines.arrow(x, lwd=lwd, col=col, lty=lty));
	invisible();
}

### Chemistry

#' @export
plot.molecule = function(x, y = NULL, lwd = 1,
		col = 1, col.ends = "red", col.id = "blue",
		adj.lim = NULL, ...) {
	if(inherits(x, "bioshape")) {
		xy1 = x[[1]][[1]];
		xy2 = if(length(x) == 1) NULL else x[-1];
		x = xy1$x;
		y = xy1$y;
	}
	# Plot:
	if(! is.null(adj.lim)) {
		if(inherits(adj.lim, "numeric")) {
			adj.lim = list(x = adj.lim[1], y = adj.lim[2]);
		}
		lim.x = adjust.range(x, adj.lim$x);
		lim.y = adjust.range(y, adj.lim$y);
		plot(x, y, type = "l", asp = 1, lwd=lwd, col=col,
			xlim = lim.x, ylim = lim.y, ...);
	} else {
		plot(x, y, type = "l", asp = 1, lwd=lwd, col=col, ...);
	}
	if( ! is.null(xy2)) lines.object.base(xy2, lwd=lwd, col=col);
	if( ! is.null(col.ends)) {
		id = c(1, length(x));
		points(x[id], y[id], col = col.ends);
	}
	if( ! is.null(col.id)) {
		text(jitter(x), jitter(y),
			labels = seq(length(x)), col = col.id);
	}
}

#' @export
lines.chemistry = function(x, lwd=NULL, col=1, ...) {
  lines.object.base(x, lwd=lwd, col=col, ...)
  invisible();
}

### Liposome
#' @export
lines.liposome = function(x, col="#48B000", col.line=1, lwd=1, lipid.border=FALSE) {
	lines.circles(x[[1]], line = lipid.border, fill=col);
	lines.circles(x[[2]], line = lipid.border, fill=col);
	# Lines: requires list: x[3]!
	lines.object.base(x[3], lwd=lwd, col=col.line);
	invisible();
}

### Circle
#' @export
lines.circle = function(x, lwd = NULL, col = NULL, ...) {
  lines.object.base(list(x), lwd=lwd, col=col, ...)
  invisible();
}


### Circles
# - collection/chain of circles;
# fill = fill colour;
# col  = colour of circle borders;
# col.line = colour of extra circle of radius R;
#' @export
lines.circles = function(x, R, fill="#B0B032", col=NULL, col.line="green", line=TRUE, ...) {
  xy = x;
  r = attr(xy, "r");
  if(is.null(r)) stop("Missing r!");
  x = xy$x; y = xy$y;
  n = length(x);
  if(length(r) < n) r = rep(r, n %/% length(r));
  #
  for(id in seq(n)) {
    if(is.null(fill)) {
      shape::filledcircle(r1=r[id], r2=0, mid = c(x[id], y[id]), ...)
    } else {
      shape::filledcircle(r1=r[id], r2=0, mid = c(x[id], y[id]), col=fill, ...);
    }
  }
  if(line) {
    R = if(is.null(R)) attr(xy, "R") else R;
    center = attr(xy, "center");
    shape::plotcircle(r=R, mid=center, lcol=col.line, col=NULL);
  }
}

### Ellipse
#' @export
lines.ellipse = function(x, lwd = NULL, col = NULL, ...) {
  lines.object.base(list(x), lwd=lwd, col=col, ...)
  invisible();
}


###############

### Other

# L = length of the line segment;
#' @export
lines.slope = function(xy, slope, L = 4, ...) {
	L = L / 2;
	# (x2 - x1)^2 = L^2 / (sl^2 + 1)
	dx = L / sqrt(slope^2 + 1);
	x2 = xy[1] + c(-dx, dx);
	dy = L / sqrt(1 + 1/slope^2);
	if(slope < 0) dy = - dy;
	y2 = xy[2] + c(-dy, dy);
	lines(x2, y2, ...);
}

######################

######################
### Plot Functions ###

### Plot Torus
#' @export
lines.torus = function(r, center = c(0,0), col = c(1,1), fill = col,
		lwd = c(1,1), lwd.fill = 1, ..., N = c(64, 32)) {
	xy1 = points.circle(N[1], r[1], center=center);
	xy2 = points.circle(N[2], r[2], center=center);
	as.xym = function(x) {
		x$x = c(x$x, x$x[1]);
		x$y = c(x$y, x$y[1]);
		cbind(x$x, x$y);
	}
	if(! is.null(fill)) {
		shape::filledshape(as.xym(xy1), as.xym(xy2), lcol = NA, col = fill, lwd = lwd.fill);
	}
	if(length(col) == 1) col = c(col, col);
	if(length(lwd) == 1) lwd = c(lwd, lwd);
	polygon(xy1, border = col[1], lwd = lwd[1], ...);
	polygon(xy2, border = col[2], lwd = lwd[2], ...);
}

