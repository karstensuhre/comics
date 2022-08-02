# comics
**The Molecular Human – A Roadmap of Molecular Interactions Linking Multiomics Networks with Disease Endpoints**

This is the GitHub page for the *comics* multiomics web-server.

*Comics* is a [shiny](https://shiny.rstudio.com) interface to the *The Molecular Human* and provides a visualization of over 34,000 associations between over 8,700 multiomics traits and disease endpoints.

The underlying networks are derived from 18 technically diverse deep molecular phenotyping (omics-)platforms analyzing urine, blood, and saliva samples from up to 374 participants of the multi-ethnic diabetes case-control study QMDiab.

The links between the multiomics traits include partial correlations between traits from the same platform (GGMs), mutual best hits for pairwise all-against-all correlations between platforms, genome-wide, epigenome-wide and transcriptome-wide multiomics and disease associations, and associations with QMDiab clinical endpoints.

Details of the *comics* server and its underlying multiomics networks can be found on the [*comics* homepage](http://www.metabolomix.com/comics/).

## Running comics via a web server 
A comics web server is available at [shinyapps.io](https://www.shinyapps.io/) following this link: http://comics.metabolomix.com.

## Running comics locally using docker
You need to download and install docker as described [here](https://www.docker.com/get-started/).

Then pull the docker image `ghcr.io/karstensuhre/comics` using the docker GUI or run the following command line:

```
docker run --detach --name comics -p 8080:3838 ghcr.io/karstensuhre/comics
```

Navigate to your browser and open the following page: http://localhost:8080.

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

Here is a screenshot of the `comics` server:

![Comics Screenshot](ComicsScreenshot1.png)
