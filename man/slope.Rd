\name{slope}
\alias{slope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Computes the slope of a line
}
\description{
	Basic utility function which computes the slope of a line.
}
\usage{
	slope(x, y)
}

\arguments{
  \item{x, y}{
	the x and y coordinates of the line segment;
}
}

\details{
The function requires at least 2 points for the x and y coordinates. If more than 2 points are provided, then the slope of each individual line segment is computed.
}
\value{
Returns the slope of the line.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{atan}} and \code{\link{atan2}} compute the angle.
}
\examples{

slope(c(1,3), c(1,3))

slope(1:3, 3:1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ ~dplot }
% \keyword{ ~kwd2 }
\concept{ ~line }
% \concept{ ~cpt2 }
