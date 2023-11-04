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


### Functions to Generate Arrow-Heads


#####################

### Simple ArrowHead: --->
# (x, y) = tip of the ArrowHead;
#' @export
arrowHeadSimple = function(x, y, slope, d=-1, dV=c(-d, d), scale=1) {
  p = if(d == 0) matrix(c(x, y), nrow=1, ncol=2)
  else shiftPoint(c(x, y), slope=slope, d = d, scale=scale);
  pV = shiftLine(p, slope=slope, d = dV, scale=scale);
  arrHead = list(
    x = c(pV[1,1], x, pV[2,1]),
    y = c(pV[1,2], y, pV[2,2]));
  return(arrHead);
}

# Diamond ArrowHead: ---<>
#' @export
arrowHeadDiamond = function(x, y, slope, d=-1, dV=c(d, -d), scale=1) {
  if(length(d) > 2) stop("Only 2 values are supported for d!");
  d1 = d[[1]];
  d2 = if(length(d) == 1) 2*d else sum(d);
  # TODO: more than 2 values for dV;
  pV = arrowHeadSimple(x, y, slope=slope, d=d1, dV=dV, scale=scale);
  p2 = shiftPoint(c(x, y), slope=slope, d = d2);
  arrHead = list(
    x = c(pV$x[1], p2[1,1], pV$x[3], x, pV$x[1]),
    y = c(pV$y[1], p2[1,2], pV$y[3], y, pV$y[1]));
  return(arrHead);
}

# Double Lined ArrowHead: --->>
# - a high-level helper function;
#' @export
arrowHeadDouble = function(x, y, slope, d=-1, dH=d, dV=c(dH, -dH), scale=1) {
  # First Arrow:
  arrHead = list(H1 = arrowHeadSimple(x, y, slope=slope, d = dH, dV=dV, scale=scale));
  # Double Arrow:
  # - firstly shift point:
  p2 = shiftPoint(c(x, y), slope=slope, d=d, scale=scale);
  arrHead2 = list(H2 = arrowHeadSimple(p2[1], p2[2], slope=slope, d = dH, dV=dV, scale=scale));
  arrHead  = c(arrHead, arrHead2);
  return(arrHead);
}

### Inverted ArrowHead: ---<
# (x, y) = potential tip of the ArrowHead (when d <= 0);
#' @export
arrowHeadInverted = function(x, y, slope, d=-1, dV=c(-d, d), isTip = TRUE, scale=1) {
  pV = c(x, y);
  p  = if(d == 0) matrix(pV, nrow=1, ncol=2)
  else shiftPoint(pV, slope=slope, d=d, scale=scale);
  #
  if(isTip) {
    pA = shiftLine(pV, slope=slope, d=dV, scale=scale);
    pV = p;
  } else {
    pA = shiftLine(p[1,], slope=slope, d=dV, scale=scale);
  }
  arrHead = list(
    x = c(pA[1,1], pV[1], pA[2,1]),
    y = c(pA[1,2], pV[2], pA[2,2]));
  return(arrHead);
}

# N-Lined ArrowHead: --->>...> (n times)
#' @export
arrowHeadN = function(x, y, slope, n=1, d = 0.5, dH = - d, dV=c(dH, -dH), scale=1) {
  # Shift point along line:
  arrHead = list(arrowHeadSimple(x, y, slope=slope, d = dH, dV = dV, scale=scale));
  if(n == 1) return(arrHead);
  # Double Arrow
  for(id in seq(n-1)) {
    p = shiftPoint(c(x, y), slope=slope, d = - id * d, scale=scale);
    arrowhead = list(arrowHeadSimple(p[1], p[2], slope=slope, d = dH, dV = dV, scale=scale));
    arrHead  = c(arrHead, arrowhead);
  }
  return(arrHead);
}

### Double Lined Inverted ArrowHead: ---<<
# dH = abs(d) ensures always inverted!
# dH = width, dV = height
# d  = dist between "<<";
#' @export
arrowHeadDoubleInverted = function(x, y, slope, d=-1, dH=0.5, dV=c(-dH, dH), scale=1) {
  # Shift point along line:
  dH = abs(dH);
  dH = if(d <= 0) - dH else dH;
  # Head: 2nd "<" of "--<<"
  pV = shiftPoint(c(x, y), slope=slope, d = dH, scale=scale);
  arrHead2 = list(arrowHeadSimple(pV[1], pV[2], slope=slope, d = - dH, dV=dV, scale=scale));
  midpoint = list(pV);
  # Head: 1st "<" of "<<"
  pV = shiftPoint(c(x, y), slope=slope, d = d + dH, scale=scale);
  arrHead  = list(arrowHeadSimple(pV[1], pV[2], slope=slope, d = - dH, dV=dV, scale=scale));
  arrHead  = c(H1 = arrHead, H2 = arrHead2);
  midpoint = c(M1 = list(pV), M2 = midpoint);
  attr(arrHead, "Mid") = midpoint;
  return(arrHead);
}

# T ArrowHead: ---|
#' @export
arrowHeadT = function(x, y, slope, d=-1, dV=c(d, -d), scale=1) {
  p  = cbind(x, y);
  if(length(dV) == 1) dV = c(0, dV);
  pV = shiftLine(p, slope=slope, d = dV, scale=scale);
  arrHead = list(
    x = c(pV[1,1], pV[2,1]),
    y = c(pV[1,2], pV[2,2]));
  return(arrHead);
}

### Measurement ArrowHead: --->|
#' @export
arrowHeadMeasure = function(x, y, slope, d=-1, dV=c(d, -d), dT=dV, scale=1) {
  arrHead = arrowHeadSimple(x, y, slope=slope, d=d, dV=dV, scale=scale);
  arrHead = list(arrHead, arrowHeadT(x, y, slope=slope, dV=dT, scale=scale));
  return(arrHead);
}

# X ArrowHead: ---X
#' @export
arrowHeadX = function(x, y, slope, d=-1, dV=c(d, -d), scale=1) {
  if(length(d) > 2) stop("Only 2 values are supported for d!");
  d2 = if(length(d) == 1) 2*d else sum(d);
  # TODO: 1 value & more than 2 values for dV;
  pB1 = c(x[1], y[1]);
  pB2 = shiftPoint(pB1, d=d2, slope=slope);
  p1 = shiftLine(pB1, d=dV, slope=slope, scale=scale);
  p2 = shiftLine(pB2, d=dV, slope=slope, scale=scale);
  if(length(d) == 1) {
    midpointX = (p1$x[2]+p2$x[1])/2;
    midpointY = (p1$y[2]+p2$y[1])/2;
    arrHead = data.frame(
      x = c(p1$x[2], p2$x[1], p2$x[2], p1$x[1]),
      y = c(p1$y[2], p2$y[1], p2$y[2], p1$y[1]),
      id = c(1,1, 2,2));
  } else {
    pM = shiftPoint(c(x[1], y[1]), d=d[1], slope=slope);
    midpointX = pM[1];
    midpointY = pM[2];
    arrHead = data.frame(
      x = c(p1$x[1], midpointX, p1$x[2], p2$x[1], midpointX, p2$x[2]),
      y = c(p1$y[1], midpointY, p1$y[2], p2$y[1], midpointY, p2$y[2]),
      id = c(1,1,1, 2,2,2));
  }
  attr(arrHead, "Mid") = c(midpointX, midpointY);
  return(arrHead);
}

### Square ArrowHead: |_|
#' @export
arrowHeadSquare = function(x, y, slope, d=-1, dV=c(d, -d), scale=1) {
  if(length(d) > 1) stop("Only 1 value is supported for d!");
  # TODO: more than 2 values for dV;
  pB1 = c(x[1], y[1]);
  pB2 = shiftPoint(c(x[1], y[1]), d=d, slope=slope, scale=scale);
  p1 = shiftLine(pB1, d=dV, slope=slope, scale=scale);
  p2 = shiftLine(pB2, d=dV, slope=slope, scale=scale);
  arrHead = list(
    x = c(p1$x[2], p1$x[1], p2$x[1],  p2$x[2], p1$x[2]),
    y = c(p1$y[2], p1$y[1], p2$y[1],  p2$y[2], p1$y[2]));
  midpoint = pB2;
  attr(arrHead, "Mid") = midpoint;
  return(arrHead);
}

# Circle ArrowHead ---O
#' @export
arrowHeadCircle = function(x, y, slope, r=0.5, scale=1) {
  center = shiftPoint(c(x, y), slope = slope, d = -r, scale=scale)
  startP = shiftPoint(c(x, y), slope = slope, d = -2*r, scale=scale)
  lst = list(r=r, center=center);
  attr(lst, "class") = c("circle", class(lst));
  lst = list(lst);
  attr(lst, "start") = list(startP, center);
  return(lst)
}

# Triangle ArrowHead: ---|>
#' @export
arrowHeadTriangle = function(x, y, slope, d=-1, dV=c(-d, d), scale=1) {
  p = if(d == 0) matrix(c(x, y), nrow=1, ncol=2)
  else shiftPoint(c(x, y), slope=slope, d = d, scale=scale);
  pV = shiftLine(p, slope=slope, d = dV, scale=scale);
  arrHead = list(
    x = c(pV[1,1], x, pV[2,1], pV[1,1]),
    y = c(pV[1,2], y, pV[2,2], pV[1,2]));
  attr(arrHead, "Mid") = p;
  return(arrHead);
}
