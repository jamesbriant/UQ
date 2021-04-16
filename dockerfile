FROM rocker/rstudio:3.6.3

# RUN install2.r mvtnorm deSolve
# RUN install2.r data.table plotly RSpectra pracma fields

RUN  echo 'install.packages(c("data.table"), \
repos="http://cran.us.r-project.org", \
dependencies=TRUE)' > /tmp/packages.R \
  && Rscript /tmp/packages.R
  
RUN install2.r plotly RSpectra pracma fields