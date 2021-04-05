FROM rocker/rstudio:3.6.3

# RUN install2.r mvtnorm deSolve
RUN install2.r data.table plotly RSpectra