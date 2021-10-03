#' A toy data set with a numbat shape hidden among noise variables
#'
#' There are 7 variables (x1-x7) and 2,100 observations.
#' Variables 4 and 7 have the numbat. The rest are
#' noise. Group A has the numbat, and group B is all noise.
#'
#' @docType data
#' @name numbat
NULL

#' Data from Anscombe's famous example in tidy format
#'
#' All variables and pairs of variables have same
#' summary statistics but are very different data,
#' as can be seen by visualisation.
#'
#' @format A tibble with 44 observations and 3 variables
#' \describe{
#'   \item{set}{label of the data set, each set has 11 observations}
#'   \item{x}{variable for horizontal axis}
#'   \item{y}{variable for vertical axis}
#' }
#' @docType data
#' @name anscombe_tidy
NULL

#' datasaurus_dozen data
#'
#' From the datasauRus package. A modern update of Anscombe.
#' All plots have same x and y mean, variance and correlation,
#' but look different visually.
#'
#' @format A tibble with 1,846 observations and 3 variables
#' \describe{
#'   \item{dataset}{label of data set}
#'   \item{x}{variable for horizontal axis}
#'   \item{y}{variable for vertical axis}
#' }
#'
#' @docType data
#' @name datasaurus_dozen
NULL

#' Data from Anscombe's famous example in tidy wide format
#'
#' All variables and pairs of variables have same
#' summary statistics but are very different data,
#' as can be seen by visualisation.
#'
#' @format A tibble with 142 observations and 26 variables
#' \describe{
#'   \item{away_x, away_y}{x and y variables for away data}
#'   \item{bullseye_x, bullseye_y}{x and y variables for bullseye data}
#'   \item{circle_x, circle_y}{x and y variables for circle data}
#'   \item{dino_x, dino_y}{x and y variables for dino data}
#'   \item{dots_x, dots_y}{x and y variables for dots data}
#'   \item{h_lines_x, h_lines_y}{x and y variables for h_lines data}
#'   \item{high_lines_x, high_lines_y}{x and y variables for high_lines data}
#'   \item{slant_down_x, slant_down_y}{x and y variables for slant_down data}
#'   \item{slant_up_x, slant_up_y}{x and y variables for slant_up data}
#'   \item{star_x, star_y}{x and y variables for star data}
#'   \item{v_lines_x, v_lines_y}{x and y variables for v_lines data}
#'   \item{wide_lines_x, wide_lines_y}{x and y variables for wide_lines data}
#'   \item{star_x, star_y}{x and y variables for star data}
#'   \item{x_shape_x, x_shape_y}{x and y variables for x_shape data}
#' }
#'
#' @docType data
#' @name datasaurus_dozen_wide
#' @rdname datasaurus_dozen
NULL

#' Simulated data with special features
#'
#' Simulated data with common features that might
#' be seen in 2D data. Variable are feature, x, y.
#'
#' @format A tibble with 1,013 observations and 3 variables,
#' and 15 different patterns
#' \describe{
#'   \item{feature}{label of data set}
#'   \item{x}{variable for horizontal axis}
#'   \item{y}{variable for vertical axis}
#' }
#'
#' @docType data
#' @name features
NULL

#' Parkinsons data from UCI machine learning archive
#'
#' Biomedical voice measurements from 31 people,
#' 23 with Parkinson's disease (PD). Each column
#' in the table is a particular voice measure, and
#' each row corresponds one of 195 voice recording
#' from these individuals ("name" column). The main
#' aim of the data is to discriminate healthy people
#' from those with PD, according to "status" column
#' which is set to 0 for healthy and 1 for PD.
#'
#' The data is available at [The UCI Machine Learning Repository](https://archive.ics.uci.edu/ml/datasets/Parkinsons)
#' in ASCII CSV format. The rows of the CSV file contain
#' an instance corresponding to one voice recording.
#' There are around six recordings per patient, the name
#' of the patient is identified in the first column.
#'
#' The data are originally analysed in:
#' Max A. Little, Patrick E. McSharry, Eric J. Hunter, Lorraine O. Ramig (2008),
#' 'Suitability of dysphonia measurements for telemonitoring of Parkinson's disease',
#' IEEE Transactions on Biomedical Engineering.
#'
#' @format A tibble with 1,013 observations and 3 variables
#' \describe{
#'   \item{name}{ASCII subject name and recording number}
#'   \item{`MDVP:Fo(Hz)`}{Average vocal fundamental frequency}
#'   \item{`MDVP:Fhi(Hz)`}{Maximum vocal fundamental frequency}
#'   \item{`MDVP:Flo(Hz)`}{Minimum vocal fundamental frequency}
#'   \item{`MDVP:Jitter`,`MDVP:Jitter(Abs)`,`MDVP:RAP`,`MDVP:PPQ`,`Jitter:DDP`}{Several measures of variation in fundamental frequency}
#'   \item{`MDVP:Shimmer`,`MDVP:Shimmer(dB)`,`Shimmer:APQ3`,`Shimmer:APQ5`,`MDVP:APQ`,`Shimmer:DDA`}{Several measures of variation in amplitude}
#'   \item{`NHR`,`HNR`}{Two measures of ratio of noise to tonal components in the voice}
#'   \item{`status`}{Health status of the subject (one) - Parkinson's, (zero) - healthy}
#'   \item{`RPDE`,`D2`}{Two nonlinear dynamical complexity measures}
#'   \item{`DFA`}{Signal fractal scaling exponent}
#'   \item{`spread1`,`spread2`,`PPE`}{Three nonlinear measures of fundamental frequency variation}
#' }
#'
#' @docType data
#' @name pk
NULL
