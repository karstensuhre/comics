set -x

# copy files to assets
rm -r -f assets
mkdir assets
mkdir assets/www
cp ../multiomics-docker/generate_network.R assets/app.R
cp ../multiomics-docker/DNA_RNA_METH_supptables_v5.xlsx assets
cp ../multiomics-docker/www/about.jpg assets/www
cp ../multiomics-docker/www/legend_grey.jpg assets/www
cp ../multiomics-docker/howto.html assets
cp ../multiomics-docker/usecases.html assets

# create the Dockerfile
cat > Dockerfile <<'EOF'
FROM rocker/shiny-verse:4.2.1

# copy files to the server
COPY ./assets/app.R /srv/shiny-server/
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
EOF

# use the new build version
DOCKER_BUILDKIT=1

# build the image
docker.exe build --no-cache -f Dockerfile -t ghcr.io/karstensuhre/comics . || exit

# stop and remove the old container
docker.exe stop comics
docker.exe container rm comics

# start the image
docker.exe run --detach --name comics -p 8080:3838 ghcr.io/karstensuhre/comics

# ssh into the image
docker.exe exec -it comics /bin/bash


