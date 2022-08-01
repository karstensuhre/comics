# comics
The Molecular Human – A Roadmap of Molecular Interactions Linking Multiomics Networks with Disease Endpoints.

## Running comics via a web server 
A comics web server is available at [shinyapps.io](https://www.shinyapps.io/) following this link: http://comics.metabolomix.com

## Running comics using docker
You need to download and install docker as described [here](https://www.docker.com/get-started/).

Then use the docker image `ghcr.io/karstensuhre/comics` using the GUI or run the following command line:

```
docker run --detach --name comics -p 8080:3838 ghcr.io/karstensuhre/comics
```

Then navigate to your browser and open the following page: http://localhost:8080

