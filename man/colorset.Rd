\name{colorset}
\alias{colorset}
\title{retrieve color setting for logo}
\description{
	retrieve color setting for logo
}
\usage{
colorset(alphabet="DNA", colorScheme='auto')
}
\arguments{
\item{alphabet}{character, 'DNA', 'RNA' or 'AA'}
\item{colorScheme}{'auto', 'charge', 'chemistry', 'classic' or 'hydrophobicity' for AA, 'auto', 'basepairing', or 'blindnessSafe' for DNA ro RNA}
}
\value{
	A character vector of color scheme
}
\examples{
col <- colorset("AA", "hydrophobicity")
}