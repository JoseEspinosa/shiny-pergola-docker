############################################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Jan 2017                                        ###
############################################################################################
### Shiny app to show pergola data                                                       ###
### ui.R                                                                                 ###
############################################################################################
### TODO                                                                                 ###
### Change color scheme with a color blind friendly scheme, some ideas here:             ###
### http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/#a-colorblind-friendly-palette ###
############################################################################################

# Running the app
# runApp ("/Users/jespinosa/git/shinyPergola/apps/minPhenBrowser")

# library(shiny)
# library("Gviz")
# library('GenomicRanges')
# library ("rtracklayer")

# fluidPage(
#   # Application title
#   titlePanel("Word Cloud"),
#   
#   sidebarLayout(
#     # Sidebar with a slider and selection inputs
#     sidebarPanel(
#       selectInput("selection", "Choose a book:",
#                   choices = books),
#       actionButton("update", "Change"),
#       hr(),
#       sliderInput("freq",
#                   "Minimum Frequency:",
#                   min = 1,  max = 50, value = 15),
#       sliderInput("max",
#                   "Maximum Number of Words:",
#                   min = 1,  max = 300,  value = 100)
#     ),
#     
#     # Show Word Cloud
#     mainPanel(
#       plotOutput("plot")
#     )
#   )
# )


shinyUI(
  fluidPage(
#     titlePanel("Behavioral browser")# ,
    headerPanel(#HTML("<strong>Behavioral browser</strong>")
      HTML('Behavioral browser
           <a href="http://cbcrg.github.io/pergola/" target="_blank"><img align="right" alt="Pergola logo" 
           src="https://cloud.githubusercontent.com/assets/6224346/12887167/dcf80b24-ce72-11e5-8389-90122fd6c84e.png" /></a>')
#                 tags$head(tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
#                           tags$style(type="text/css", "select { max-width: 200px; }"),
#                           tags$style(type="text/css", "textarea { max-width: 185px; }"),
#                           tags$style(type="text/css", ".jslider { max-width: 200px; }"),
#                           tags$style(type='text/css', ".well { max-width: 330px; }"), # left menu
#                           tags$style(type='text/css', ".span4 { max-width: 330px; }")) 
    ),
    
    sidebarPanel(
      # Load info such as experimental info
      conditionalPanel(condition="input.tabs_p=='Browser'",
#                        uiOutput("windowsize"), # server.R
                       uiOutput("bedGraphRange"),
                       uiOutput("dataInterval"),
                       uiOutput("idSelect"),                       
#                        uiOutput("genomicPositionSelect"),
                       uiOutput("groups"),
                       uiOutput("groups_plot"),
                       uiOutput("boxplot")
      ),
            conditionalPanel(condition="input.tabs_p=='About'",
                             h4("Introduction") 
            )#,
      #       conditionalPanel(condition="input.tabs_p=='Plots'",
      #                        h4("Plots") 
      #       )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Browser",
                 fluidRow(                  
                  column(12,
                           plotOutput("plotbed", height=1200),
                           plotOutput("legend_track", height=30),
                           plotOutput("envInfo", height=20))#,
#                           textOutput("text1")
                        ),
                  column(3, downloadButton("all_plot_tiff", "Download snapshot"))
                ),        
        #         tabPanel("Plots",
        #                  
        #                  fluidRow(
        #                    column(3, downloadButton("barPlotValueTiff", "Download tiff")),
        #                    column(9, plotOutput("barPlotValue"))
        #                  ),
        #                  fluidRow(
        #                    column(3, downloadButton("barPlotDurationTiff", "Download tiff")),
        #                    column(9, plotOutput("barPlotDuration"))
        #                  ),
        #                  #output$barPlotNumberTiff 
        #                  fluidRow(
        #                    column(3, downloadButton("barPlotNumberTiff", "Download tiff")),
        #                    column(9, plotOutput("barPlotN"))
        #                  ),
        #                  fluidRow(
        #                    column(3, downloadButton("barPlotRateTiff", "Download tiff")),
        #                    column(9, plotOutput("barPlotRate"))
        #                  )
        #         ),        
        #         tabPanel("Upload data",
        #                  HTML('<p>Additional data</p>'),  
        #                  tableOutput('bed'),     
        #                  HTML('<p>Phases data</p>'),  
        #                  tableOutput('fileEnv'),
        #                  HTML('<p>End data</p>')),    
        tabPanel("About",
                 HTML('<p>Pergola</p>')),
        
        id="tabs_p"
        
      )
    )
  )
)