# shiny-pergola

Dockerized Shiny Pergola app

## Introduction

Data processed using [Pergola](http://cbcrg.github.io/pergola/) can be
interactively explore using a [Shiny](https://shiny.rstudio.com/) based app.

## Getting started

### Dependencies

shiny-pergola only needs Docker installed in your machine 

    Docker 1.0

### Setup

Download shiny-pergola Docker container

```
docker pull pergola/shiny-pergola:0.1
```

Run dockerized version of shiny-pergola with your data produced using Pergola:

```
docker run --rm -p 80:80 -v /your/local/path/to/pergola/folder:/pergola_data pergola/shiny-pergola:0.1 &
```









