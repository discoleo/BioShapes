#######################################
#
# BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis (2022-2023)
# Candidate: Adrian Cotoc
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
# GitHub: https://github.com/Adi131313/BioShapes
#
# 2. Bachelor Thesis: Darian Voda (2021-2022)


#### Cell-like Objects ####


### Smooth muscle cell or Fibroblast
#' @export
cell.SmoothMuscle = function(x, y, r = 1, r.nc = r/2, slope = NULL, t.nc = c(0.55, 0.48),
		lwd = 1, col = NULL, fill = NULL, col.nc = 1, fill.nc = NULL, N=128, phi = c(0, pi)) {
	if(is.null(slope)) slope = slope(x, y);
	d = sqrt((x[1] - x[2])^2 + (y[1] - y[2])^2);
	dx = d / N;
	pL = seq(0, d, by=dx);
	pp = shiftPoint(c(x[1], y[1]), d=pL, slope=slope);
	pL = pL + x[1];
	px = pp[,1] - x[1]; px = px * pi / max(abs(px));
	# Margin 1:
	pS = r * sin(px + phi[1]) + pp[,2];
	lst = list(x = pL, y = pS);
	# Margin 2:
	pS = r * sin(px + phi[2]) + pp[,2];
	lst$x = c(lst$x, rev(pL));
	lst$y = c(lst$y, rev(pS));
	#
	lst$lwd = lwd;
	if( ! is.null(col)) lst$col = col;
	if( ! is.null(fill)) lst$fill = fill;
	class(lst) = c("polygon", "list");
	# Nucleus:
	if(r.nc != 0) {
		# TODO: ellipse;
		mid.x = (1-t.nc[1])*x[1] + t.nc[1]*x[2];
		mid.y = (1-t.nc[2])*y[1] + t.nc[2]*y[2];
		center = c(mid.x, mid.y);
		nc = list(r = r.nc, center = center, col = col.nc);
		if( ! is.null(fill.nc)) nc$fill = fill.nc;
		class(nc) = c("circle", "list");
		lst = list(Cell = lst, N = nc);
	} else lst = list(lst);
	lst = as.bioshape(lst);
	return(lst)
}


### Brush-Border Cells
# p1 = Base of Cell (Point 1);
# w, h = width, height;
# n = number of half-cycles of sine-wave;
#' @export
as.options.BrushBorder = function(...) {
	opt = list(...);
	len = length(opt);
	if(len == 1 && is.list(opt[[1]])) opt = opt[[1]];
	if(is.null(opt$n)) opt$n = 6.5;
	if(is.null(opt$A)) opt$A = 0.5;
	if(is.null(opt$N)) opt$N = 128;
	if(is.null(opt$phi)) opt$phi = 0;
	return(opt);
}
#' @export
cell.BrushBorder = function(p1, w, h, slope=0, r.nc = ~ 1/5, t.nc = c(1/2, 7/20),
		lwd = 1, lwd.nc = lwd, col = NULL, fill = NULL, col.nc = 1, fill.nc = NULL,
		brush.options = list(n=6.5, A=0.5, N=128, phi=0)) {
	# Cell:
	p11 = p1;
	p12 = shiftPoint( p1, d=w, slope=slope);
	p21 = shift.ortho(p1, d=h, slope=slope);
	p21 = unlist(p21);
	p22 = shiftPoint(p21, d=w, slope=slope);
	# Brush-Border:
	opt   = as.options.BrushBorder(brush.options);
	brush = helix(p21, p22, n = opt$n, A = opt$A, phi = opt$phi, N = opt$N);
	brush = brush[[1]];
	brush$x = c(brush$x, p22[1], p12[1], p11[1], p21[1]);
	brush$y = c(brush$y, p22[2], p12[2], p11[2], p21[2]);
	brush$lwd = lwd;
	if( ! is.null(col)) brush$col = col;
	if( ! is.null(fill)) brush$fill = fill;
	class(brush) = c("polygon", "list");
	lst = list(Cell = brush);
	# Nucleus:
	isFormula = inherits(r.nc, "formula");
	if(isFormula || r.nc != 0) {
		t.nc[2] = 1 - t.nc[2];
		center = center.p4(p11, p12, p21, p22, t = t.nc);
		if(isFormula) {
			d = sqrt((p11[1] - p12[1])^2 + (p11[2] - p12[2])^2);
			r = eval(r.nc[[2]]) * d;
		} else r = r.nc;
		nc = list(r = r, center = center, lwd = lwd.nc, col = col.nc);
		if( ! is.null(fill.nc)) nc$fill = fill.nc;
		class(nc) = c("circle", "list");
		lst$N = nc;
	}
  return(as.bioshape(lst));
}

### Muscle tissue ###
#' @export
muscle = function(scale.x = 1.5, scale.R = c(1.5, 1.5),
                  x = c(-2, 2), y = c(1, 1), dy = 0.4, dx = 2, n = 6, fill = "red"){

  lst = list()
  for(iy in seq(0,n)){
    if(iy%%2 == 0){
      for(idx in seq(-2,2,2)){
        lens = lens(x = x + idx*dx, y = y + iy*dy, scale=scale.R, fill = fill)
        lst = c(lst, lens)
      }}
    else{
      for(idx in seq(-3,3,2)){
        lens = lens(x = x + idx*dx, y = y + iy*dy, scale=scale.R, fill = fill)
        lst = c(lst, lens)
      }
    }}
  lst = as.bioshape(lst)
  return (invisible(lst))
}

### Red blood cells
#' @export
draw_blood_cell = function(center = c(0, 0),
                           radius = 1, lwd = 10, col = "#901000", fill = "red"){
  #lim = c(-radius, radius) * 2;
  #plot.base(xlim = lim, ylim = lim)
  #circle <- seq(0, 2 * pi, length.out = 80)
  #x <- radius * sin(circle)
  #y <- radius * cos(circle)
  #polygon(x, y, col = "red", border = "#901000", lwd = lwd)
  lst = list(center = center, r = radius,
             fill = fill, col = col, lwd = lwd);
  class(lst) = c("circle", "list");
  lst = list(lst);
  class(lst) = c("bioshape", "list");
  return(lst)

}
