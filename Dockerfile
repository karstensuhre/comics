FROM rocker/shiny-verse:4.2.1

# copy files to the server
COPY ./assets/generate_network.R /srv/shiny-server/app.R
COPY ./assets/DNA_RNA_METH_supptables_v5.xlsx /srv/shiny-server/
COPY ./assets/www/about.jpg /srv/shiny-server/www/
COPY ./assets/www/legend_grey.jpg /srv/shiny-server/www/
COPY ./assets/howto.html /srv/shiny-server/
COPY ./assets/usecases.html /srv/shiny-server/

# install ubuntu software
RUN apt-get update && apt-get install -y --no-install-recommends xdg-utils vim

# install R packages
RUN install2.r --error --skipinstalled \
    visNetwork \
    shinydashboard \
    DT

USER shiny
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
