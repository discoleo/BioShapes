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


### Examples: Arrows & Labels


### Arrows: Summary
#' @export
example.arrows = function(dx = c(0, 0), lwd = 2, lty = 1, d.lines = 0, join = 0,
		fill = c("red", "#FFA090", "#3296F2", "#E648D0"),
		new.plot = TRUE) {
  if(length(join) == 1) join = rep(join, 3);
  cex = 0.75; cexSmall = 0.70;
  
  ### Plot
  if(new.plot) plot.base();

  ### Row 1:
  y = c(7, 9); yt = 6.7;
  abline(h = y[2], col="green")
  text(4, y[2] + 1,
       "Arrows with different ArrowHeads")

  # Simple ArrowHead
  x = c(-1, -1) + dx;
  d = -0.5;
  d.head = c(-0.5,0.5)
  a1 = arrowSimple(x, y, lwd=lwd, lty=lty, d=d, d.head=d.head, d.lines=d.lines);
  text(-1, yt,
       "Simple", cex=cex);

  # Double ArrowHead
  # Note: use join = 1 to keep the lines non-joined;
  # - join = 0: use default = joined;
  x = c(1, 1) + dx;
  d = -0.3;
  d.head=-0.5
  a2 = arrowDouble(x, y, lwd=lwd, lty=lty, d=d, d.head=d.head, d.lines=d.lines, join=join[[1]]);
  text(1, yt,
       "Double", cex=cex);

  # Inverted ArrowHead
  x = c(3, 3) + dx;
  d = -0.5;
  d.head=c(-0.5,0.5)
  a3 = arrowInverted(x, y, lwd=lwd, lty=lty, d=d, d.head=d.head, d.lines=d.lines);
  text(3, yt,
       "Inverted", cex=cex);

  # Double Lined Inverted ArrowHead
  x = c(5, 5) + dx;
  arrowDoubleInverted(x, y, lwd=lwd, lty=lty, d=-0.3, d.lines=d.lines, join=join[[2]]);
  text(5, yt, adj = c(0.5, 0.75),
       "Double-Lined\nInverted", cex = cexSmall);

  # T Shape ArrowHead
  x = c(7, 7) + dx;
  arrowT(x, y, lwd=lwd, lty=lty, d=-0.75, d.lines=d.lines);
  text(7, yt,
       "T Shape", cex=cex);

  # Measurement ArrowHead
  x = c(9.5, 9.5) + dx;
  arrowMeasure(x, y, lwd=lwd, lty=lty, d=-0.5, d.lines=d.lines);
  text(9.5, yt,
       "Measurement", cex = cexSmall);


  ### Row 2
  y = c(3, 5); yt = 2.7;
  abline(h = y[2], col="green")

  # X Shape ArrowHead
  x = c(-1, -1) + dx;
  arrowX(x, y, lwd=lwd, lty=lty, d=0.5, d.lines=d.lines);
  text(-1, yt,
       "X Shape", cex=cex);

  # Square Shape ArrowHead
  x = c(1, 1) + dx;
  arrowSquare(x, y, lwd=lwd, lty=lty, d=-0.5, d.lines=d.lines);
  text(1, yt,
       "Square Shape", cex = cexSmall);

  # Rectangular Flag: (2*d) x d
  # Note: use c(0, d) for real square;
  x = c(3, 3) + dx;
  arrowSquare(x, y, lwd=lwd, lty=lty, d=-0.5, d.head=c(0, 2*d), d.lines=d.lines);
  text(3.5, yt, adj = c(0.5, 0.75),
       "Rectangular\nFlag", cex = cexSmall);

  # Diamond ArrowHead
  x = c(6, 6) + dx;
  d.head = c(-0.5, 0.5);
  d = -0.5;
  arrowDiamond(x, y, lwd=lwd, lty=lty, d=d, d.head=d.head, d.lines=d.lines, join=0);
  text(6, yt,
       "Diamond", cex=cex);

  # Multiple-Lined ArrowHead
  n = 3; d = 0.5;
  x = c(8, 8) + dx;
  arrowN(x, y, n=n, lwd=lwd, lty=lty, d=d, d.lines=d.lines, join = join[[3]]);
  text(8, yt,
       "Multiple-Lined", cex = cexSmall);


  ### Row 3
  y = c(-1, 1); yt = -1.3;
  x = c(-1, -1) + dx;
  abline(h = y[2], col="green")

  # Rectangle ArrowHead: Solid
  arrowSolidSquare(x, y, lwd=lwd, lty=lty, d=-0.5, d.lines=d.lines,
		col = "darkred", fill = fill[[1]]);
  text(x[1], yt, adj = c(0.5, 0.75),
       "Solid\nRectangle", cex=cex);

  # Diamond ArrowHead: Solid
  x = x + 2;
  d.head = c(-0.5, 0.5);
  d = -0.5;
  arrowDiamond(x, y, lwd=lwd, lty=lty, d=d, d.head=d.head, d.lines=d.lines,
		fill=fill[[2]], join=0);
  text(x[1], yt, adj = c(0.5, 0.75),
       "Solid\nDiamond", cex=cex);

  # Triangle ArrowHead
  x = x + 2;
  d = -0.5;
  a1 = arrowTriangle(x, y, lwd=lwd, lty=lty, d=d, d.lines=d.lines);
  text(x[1], yt,
       "Triangle", cex=cex);

  # Solid Circle ArrowHead
  x = x + 2;
  arrowCircle(x, y, lwd=lwd, lty=lty, r=0.5, d.lines=d.lines, fill = fill[[3]]);
  text(x[1], yt, adj = c(0.5, 0.75),
       "Solid\nCircle", cex=cex);

  # Simple Circle ArrowHead
  x = x + 2;
  arrowCircle(x, y, lwd=lwd, lty=lty, r=0.5, d.lines=d.lines);
  text(x[1], yt, adj = c(0.5, 0.75),
       "Simple\nCircle", cex=cex);

  # Half Circle ArrowHead
  x = x + 2;
  arrowHalfCircle(x, y, lwd=lwd, lty=lty, r=0.5, d.lines=d.lines, fill = fill[[4]]);
  text(x[1], yt, adj = c(0.5, 0.75),
       "Half\nCircle", cex=cex);
}


### Double Halves
#' @export
test.arrow.Half = function(lwd = 2, col = c(2,3,4)) {
	plot.base()
	arrowDHalf(c(2, 6), c(3, 3), lwd=lwd, col = col[1], plot = TRUE)
	arrowDHalf(c(2, 6), c(4, 7), lwd=lwd, col = col[2], plot = TRUE)
	arrowDHalf(c(1, 1), c(3, 7), lwd=lwd, col = col[3], plot = TRUE)
}


########################
########################

### Boxes & Labels

### Basic Test
# scale = 1 / aspect ratio, used in 1 test;
#' @export
test.box.cap = function(scale = 2) {
	par.old = par(mfrow = c(2,3))
	
	plot.base()
	tmp = box.cap(c(1,7), c(2,7), fill = "red", col="#FFA0A0")
	lines(tmp)
	tmp = box.cap(c(7,1), c(4,-1), fill = "red", col="#FFA0A0")
	lines(tmp)
	
	plot.base()
	tmp = box.cap(c(1,7), c(4,4), fill = "red", col="#FFA0A0")
	lines(tmp)
	tmp = box.cap(c(7,1), c(2,2), fill = "red", col="#FFA0A0")
	lines(tmp)
	lines(arrowSimple(c(2,6), c(6,6)))
	lines(arrowSimple(c(6,2), c(0,0)))
	
	plot.base(asp=1/scale)
	tmp = box.cap(c(1,7), c(2,7), fill = "red", col="#FFA0A0", scale=scale)
	lines(tmp)
	tmp = box.cap(c(7,1), c(3,-2), fill = "red", col="#FFA0A0", scale=scale)
	lines(tmp)
	text(4, 16, paste0("Scale = ", scale, " : 1"), cex = 1.75)
	
	### Row 2:
	plot.base()
	tmp = box.cap(c(1,7), c(2,7))
	lines(tmp)
	tmp = box.cap(c(1,7), c(0,5))
	lines(tmp)
	
	plot.base()
	tmp = box.cap(c(1,7), c(7,2))
	lines(tmp)
	tmp = box.cap(c(7,1), c(0,5))
	lines(tmp)
	
	plot.base()
	tmp = box.cap(c(2,2), c(2, 7))
	lines(tmp)
	tmp = box.cap(c(5,5), c(7,2))
	lines(tmp)
	
	par(par.old)
	invisible()
}


### Elliptic Cap
# scale = 1 / aspect ratio;
#' @export
test.box.capEllipse = function(y.rel = 1/4, lwd = 1, col = "black", fill = NULL, scale = 1) {
	boxf = function(x, y) {
		box.capEllipse(x, y, y.rel=y.rel, lwd=lwd, col=col, fill=fill, scale=scale);
	}
	plot.base(asp = 1 / scale);
	lines(boxf(c(1, 4), c(3.5, 0)))
	lines(boxf(c(4, 1), c(3.5, 0)))
	# for visibility:
	dy0 = if(scale <= 1) 0 else 3 / scale;
	dy2 = 3*dy0; dy1 = dy2 + dy0;
	lines(boxf(c(0, 5), c(7, 7) + dy1))
	lines(boxf(c(5, 0), c(5, 5) + dy2))
	#
	lines(boxf(c(6, 6), c(0, 4) - 0.25))
	lines(boxf(c(8, 8), c(4, 0) - 0.25))
}


####################

### Pins

# Basic Test:
#' @export
test.pins = function(theta = c(2,1,2,1,1,8) * pi/6) {
	par.old = par(mfrow = c(3,2))
	
	plot.base()
	lines(pin.center(c(1,5), c(1,5), theta = theta[1], col.circle="red"))
	points(5,5)
	
	plot.base()
	lines(pin.center(c(1,5), c(1,5), theta = theta[2], col.circle="red", fill="#B280F2", lwd=2))
	points(5,5)
	
	plot.base()
	lines(pin.center(c(1,5), c(5,1), theta = theta[3], col.circle="red"))
	points(5,1)
	
	plot.base()
	lines(pin.center(c(5,1), c(5,1), theta = theta[4], col.circle="red"))
	points(1,1)
	
	plot.base()
	lines(pin.center(c(5,1), c(1,5), theta = theta[5], col.circle="red"))
	points(1,5)
	
	# TODO: theta >= pi ???
	plot.base()
	lines(pin.center(c(1,5), c(1,5), theta = theta[6], col.circle="red"))
	points(5,5)
	
	par(par.old)
}
