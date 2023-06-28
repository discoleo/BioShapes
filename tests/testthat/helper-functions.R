aspect_ratio = 1.4
aspect_ratio_max = 2.0
# aspect_ratio_max = 1.0

#### Helper Functions ####

Dsquare = function(xy, x0, y0) {
  sum((xy$x - x0)^2, (xy$y - y0)^2)
}

dist.polygon = function(xy) {
  xc = xy$x; xc = c(xc[-1], xc[1]);
  yc = xy$y; yc = c(yc[-1], yc[1]);
  dd = sum((xy$x - xc)^2, (xy$y - yc)^2);
  return(dd);
}

testArrow = function(h, d, dV=c(-d, d)){
  len = length(h$x);
  if(len == 3) {
    dd = Dsquare(h, h$x[2], h$y[2]) - 2*d^2 - sum(dV^2);
  } else if(len == 5) {
    dd = dist.polygon(h) - 4*d^2 - 2*sum(dV^2);
  } else stop("Not yet implemented!")
  stopifnot(round(dd, 8) == 0)
}

linesAid = function(..., id=c(1,3)) {
  h = list(...);
  lapply(h, linesAid1, id=id);
  invisible();
}
linesAid1 = function(h, id = c(1,3), col="green") {
  lines(h$x[id], h$y[id], col=col)
}
