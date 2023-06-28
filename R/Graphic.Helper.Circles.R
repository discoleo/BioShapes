###################
#
# Thesis
#
# Title: BioShapes
#
# BSC Candidate: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#
# Based on BSC Thesis: Darian Voda (2021-2022)
# Faculty of Mathematics and Informatics, UVT
#
# in collaboration with Syonic SRL
#
# GitHub: https://github.com/Adi131313/Licenta


### Helper Functions to Generate Circles


### Plot Objects formed from circles;
# - convenience function;
#' @export
testFilledCircle = function(xy, r=NULL, R=NULL, lim=NULL, line=TRUE,
                            col="#B0B032", col.line="green", add=FALSE, pin = FALSE, ...) {
  if(is.null(r)) {
    r = attr(xy, "r");
    if(is.null(r)) stop("Missing r!");
  } else {
    attr(xy, "r") = r;
  }
  ### New Plot
  if( ! add) {
    x = xy$x; y = xy$y;
    mid = attr(xy, "center");
    if(is.null(lim)) {
      R0  = attr(xy, "R");
      lim = R0 + r + 1;
      lim = c(-lim, lim);
    } else if(length(lim) == 1) {
      lim = c(-lim, lim);
      mid = c(0, 0); # remove center offset;
    } else {
      mid = c(0, 0); # remove center offset;
    }
    plot(x, y, xlim = lim + mid[1], ylim = lim + mid[2], asp = 1);
  }
  if(pin){
    pin = mean(par("pin")) + 0.25;
    par.old = par(pin = c(pin, pin));
    lines.circles(xy, R=R, line=line, fill=col, col.line=col.line, ...)
    par(par.old);
  }
  else {
    lines.circles(xy, R=R, line=line, fill=col, col.line=col.line, ...)
  }
}


#####################

### Helper Functions

# r = radius;
# phi = rotation (counter-clockwise);
#' @export
pointsCircle = function(n, r, center = c(0,0), phi=0) {
  x = r * cos(seq(0, n-1) * 2*pi/n + phi) + center[1];
  y = r * sin(seq(0, n-1) * 2*pi/n + phi) + center[2];
  lst = list(x=x, y=y);
  attr(lst, "R") = r;
  attr(lst, "center") = center;
  return(lst);
}


#### Tangent circles forming a large circle ####
#' @export
circlesOnCircle = function(n, r, center = c(0,0), phi=0) {
  R  = r / sin(pi/n);
  xy = pointsCircle(n, r=R, center=center, phi=phi);
  attr(xy, "R") = R;
  attr(xy, "r") = r;
  return(xy);
}

##### Outside a large circle of given radius #####
#' @export
circlesOutsideFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r / (1/sin(pi/n) - 1);
  R1 = r + r1;
  xy = pointsCircle(n, r=R1, center=center, phi=phi);
  attr(xy, "R") = R1;
  attr(xy, "r") = r1;
  return(xy);
}

#### Tangent circles ####

##### Forming large circle of given radius #####
#' @export
circlesOnFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r * sin(pi/n);
  xy = pointsCircle(n, r=r, center=center, phi=phi);
  attr(xy, "R") = r+r1; # reuse same attribute name ???
  attr(xy, "r") = r1;
  return(xy);
}

##### Inside a large circle of given radius #####
#' @export
circlesInFixedCircle = function(n, r, center = c(0,0), phi=0) {
  r1 = r / (1 + 1/sin(pi/n));
  R  = r - r1;
  xy = pointsCircle(n, r=R, center=center, phi=phi);
  attr(xy, "R") = R;
  attr(xy, "r") = r1;
  return(xy);
}

# Generates the points on an arc of circle.
#' @export
circle.arc = function(r=1, phi, center = c(0,0), N = 64) {
  id = seq(0, N) / N;
  dp = phi[2] - phi[1];
  dw = id * dp + phi[1];
  x  = r * cos(dw) + center[1];
  y  = r * sin(dw) + center[2];
  xy = list(x=x, y=y);
  return(xy);
}


