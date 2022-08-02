# comics
**The Molecular Human – A Roadmap of Molecular Interactions Linking Multiomics Networks with Disease Endpoints**

This is the GitHub page for the *comics* multiomics web-server.

*comics* is a [shiny](https://shiny.rstudio.com) interface to the *The Molecular Human* and provides a visualization of over 34,000 correslations and associations between over 8,700 multiomics traits and disease endpoints.

The underlying networks are derived from 18 technically diverse deep molecular phenotyping (omics-)platforms analyzing urine, blood, and saliva samples from up to 374 participants of the multi-ethnic diabetes case-control study QMDiab.

Details of the *comics* server and its underlying multiomics networks can be found on the [*comics* homepage](http://www.metabolomix.com/comics/).

## Running comics via a web server 
A comics web server is available at [shinyapps.io](https://www.shinyapps.io/) following this link: http://comics.metabolomix.com

## Running comics using docker
You need to download and install docker as described [here](https://www.docker.com/get-started/).

Then use the docker image `ghcr.io/karstensuhre/comics` using the GUI or run the following command line:

```
docker run --detach --name comics -p 8080:3838 ghcr.io/karstensuhre/comics
```

Then navigate to your browser and open the following page: http://localhost:8080

## Running comics locally in RStudio
All files required to run comics locally using [RStudio](https://www.rstudio.com) are in the `./assets` directory. Start-up `rstudio` and then launch `app.R`. 

Note: you need to have the following libraries installed in RStudio:
```
install.packages("readxl")
install.packages("tidyverse")
install.packages("visNetwork")
install.packages("shiny")
install.packages("shinydashboard")
install.packages("DT")
```
