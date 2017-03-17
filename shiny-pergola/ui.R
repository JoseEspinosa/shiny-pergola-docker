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
# shiny::runApp('git/shiny-pergola-docker/shiny-pergola')

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
      conditionalPanel(condition="input.tabs_p=='Browser'",
                       uiOutput("bedGraphRange_tab"),
                       uiOutput("dataInterval_tab"),
                       uiOutput("plots2show_tab"),
                       uiOutput("groups_tab"),
                       uiOutput("type_gr_plot_tab")
      ),
            conditionalPanel(condition="input.tabs_p=='About'",
                             h4("Introduction") 
            )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Browser",
                 fluidRow(                  
                  column(12,
                           plotOutput("plotbed", height=800),
                           plotOutput("legend_track", height=100),
                           plotOutput("envInfo", height=20))#,
#                            textOutput("text1")
                        ),
                  column(3, downloadButton("all_plot_tiff", "Download snapshot"))#,
                ),            
        tabPanel("About",
                 HTML('<p>Pergola</p>')),
        
        id="tabs_p"
        
      )
    )
  )
)