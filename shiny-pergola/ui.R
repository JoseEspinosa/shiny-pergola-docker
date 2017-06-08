############################################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Jan 2017                                        ###
############################################################################################
### Shiny app to show pergola data                                                       ###
### ui.R                                                                                 ###
############################################################################################
### TODO                                                                                 ###
### Create a nicer graphical interface                                                   ###
### Inspiration to avoid remaking until change is final in reactive                      ###
### https://stackoverflow.com/questions/31051133/how-do-i-make-sure-that-a-shiny-reactive-plot-only-changes-once-all-other-reacti
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
                       h4("Links"),                      
                       helpText(a("Get Shiny-Pergola source code on GitHub!", 
                                  href = "https://github.com/JoseEspinosa/shiny-pergola-docker", 
                                  target = "_blank")),
                        helpText(a("Pergola documentation", 
                                 href = "http://cbcrg.github.io/pergola/", target = "_blank"))                               
      ), width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Browser",
                 fluidRow(                  
                  column(12,
                           plotOutput("plotbed", height=800),
                           plotOutput("legend_track", height=50),
                           plotOutput("envInfo", height=20))#,
#                            textOutput("text1")
                        ),
                  column(3, downloadButton("all_plot_tiff", "Download snapshot"))#,
                ),            
        tabPanel("About",
                 HTML('<h4>Introduction</h4> 
                      <p>Shiny-Pergola is a Shiny application for the visualization of longitudinal 
                      behavioral data. Behavioral data should be first be processed using Pergola. 
                      Once data has been processed you can render it using a configuration file as
                      explained here.</p>')),
        
        id="tabs_p"
        
      )
    )
  )
)