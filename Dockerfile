
FROM pergola/shiny-pergola-base:0.1

MAINTAINER Jose Espinosa-Carrasco <espinosacarrascoj@gmail.com>

COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf

COPY /shiny-pergola /srv/shiny-server/

EXPOSE 80

COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]
