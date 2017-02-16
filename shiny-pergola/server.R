############################################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Jan 2017                                        ###
############################################################################################
### Shiny app to show pergola data using Gviz                                            ###
### server.R                                                                             ###
############################################################################################
### TODO                                                                                 ###
### Benchmark using system.time, benchmark library                                       ###
### Try to load plots at the beginning less time                                         ### 
### Generate bedgraph files like bed files and and object with all the bedGraph files    ###
### like know for the group plot                                                         ###
############################################################################################

library(Gviz)
library(GenomicRanges)
library (rtracklayer)

col_back_title="brown"
tr_sum_size=20

col_gr_1 <- "darkblue"
col_gr_2 <- "brown"
col_ctrl <- col_gr_1
col_case <- col_gr_2
cb_palette <- c("#999999", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7")

# Avoid problems if user set many groups
cb_palette <- rep (cb_palette, 10) 

# base_dir <- "/Users/jespinosa/git/shinyPergola/data"
# base_dir <- "/Users/jespinosa/git/shinyPergola/data/worm_data"
# base_dir <- "/Users/jespinosa/2017_phecomp_marta"
# base_dir <- "/Users/jespinosa/git/shinyPergola/data/HF_experiment"
base_dir <- "/pergola_data"

# data_dir <- dir(file.path(base_dir,"bed4test"))
# data_dir <- file.path(base_dir,"bed4test")
# data_dir <- file.path(base_dir, "GB_indidividual_files")
# data_dir <- file.path(base_dir, "bed4test_all")
data_dir <- file.path(base_dir, "files")

# exp_design_f <- "exp_info_test.txt"
exp_design_f <- "exp_info.txt"

b2v <- exp_info <- read.table(file.path(base_dir, exp_design_f), header = TRUE, stringsAsFactors=FALSE)

# exp_info$sample
perg_bed_files <- sapply(exp_info$sample, function(id) file.path(data_dir, paste(id, ".bed", sep="")))

# s2c <- exp_info <- read.table(file.path(base_dir, "exp_info_test.txt"), header = TRUE, stringsAsFactors=FALSE)
b2v <- dplyr::mutate(b2v, path = perg_bed_files)
# s2c

perg_bedg_files <- sapply(exp_info$sample, function(id) file.path(data_dir, paste(id, ".bedGraph", sep="")))

bg2v <- exp_info <- read.table(file.path(base_dir, exp_design_f), header = TRUE, stringsAsFactors=FALSE)
bg2v <- dplyr::mutate(bg2v, path = perg_bedg_files)

# unique(exp_info$condition)
# grps
g_min_start <- 100000000
g_max_end <- -100000000
min_v <- 0
max_v <- 0

l_gr_color <- mapply(function(x, col) list(col),
       unique(b2v$condition), cb_palette[1:length(unique(b2v$condition))])
# l_gr_color[["control"]]

bed2pergViz <- function (df, gr_df, format_f="BED") {
  grps <- as.character(gr_df[[setdiff(colnames(gr_df), 'sample')]])
  
  r <- lapply(unique(grps),
         function(g) {
           gr_samps <- grps %in% g
           gr_files <- df$path[gr_samps]
           
           lapply(gr_files, function (bed) {             
             id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", bed)
             bed_GR <- import(bed, format = format_f)
             min_start <- min(start(bed_GR))
             max_end <- max(end(bed_GR))
             
             if (format_f == "BED") {                              
               tr <- AnnotationTrack(bed_GR, name = paste ("", id, sep=""),
                                     fill=l_gr_color[[g]], 
                                     background.title = l_gr_color[[g]], col=NULL)#, fill=col_ctrl, background.title = col_ctrl)               
             }
             
             if (format_f == "bedGraph") {
               min_v <<- floor (min(bed_GR$score))
               max_v <<- ceiling (max(bed_GR$score))
#                scores <- as.vector(mcols(bed_GR))
#                tr <- DataTrack(bed_GR, name = paste ("", id, sep=""))#, fill=col_ctrl, background.title = col_ctrl)               
               tr <- bed_GR
             }

             if (g_min_start > min_start) { g_min_start <<- min_start }
             if (g_max_end < max_end) { g_max_end <<- max_end }

             return (tr) })
         })
  
  names(r) <- unique(grps)
  return (r)
}

l_gr_annotation_tr_bed <- bed2pergViz (b2v, exp_info)

# list_all <- list()
# list_all
# l_gr_annotation_tr_bed[c("case", "control")]

# for (i in 1:length(l_gr_annotation_tr_bed)){
#   list_gr <- lapply (l_gr_annotation_tr_bed[[i]], function (l, color=cb_palette[i]) { 
#     displayPars(l) <- list(fill=color, background.title = color, col=NULL) # coll null for boxes lines
#     return (l)
#   })
#   
#   list_all <- append(list_all, list_gr)  
# }

l_granges_bg <- bed2pergViz (bg2v, exp_info, "bedGraph") 

l_gr_data_tr_bg_tmp <- lapply (seq_along(l_granges_bg), function (i_group_exp) {
        lapply (l_granges_bg[[i_group_exp]],  function (granges_obj) {
                                                     d_track <- DataTrack(granges_obj,
                                                                type="heatmap", ylim = c(0, 0.5),
                                                                background.title = l_gr_color[[i_group_exp]],
                                                                gradient=c('white','blue'))
                                                     return (d_track)
                                                    })
     })

names (l_gr_data_tr_bg_tmp) <- names(l_granges_bg)
l_gr_data_tr_bg <-l_gr_data_tr_bg_tmp
l_gr_annotation_tr_bg <- l_gr_data_tr_bg
list_all_bg <- l_gr_data_tr_bg

  
# setdiff(l_gr_annotation_tr_bg[[1]][[1]], l_gr_annotation_tr_bg[[1]][[2]])
# subsetByOverlaps (l_gr_annotation_tr_bg[[1]][[1]], l_gr_annotation_tr_bg[[1]][[2]])
# list_all_bg <- list()
# group_lab <- c()
# color_by_tr <- c()
# for (i in 1:length(l_gr_annotation_tr_bg)){
#   group_lab <- append(group_lab, rep (names(l_gr_annotation_tr_bg)[i], length(l_gr_annotation_tr_bg[[i]])))
#   color_by_tr <- append(color_by_tr, cb_palette[i], length(l_gr_annotation_tr_bg[[i]]))
#   
#   for (j in 1:length(l_gr_annotation_tr_bg[[i]])){
#     GR <- l_gr_annotation_tr_bg[[i]][[j]]
# 
#     id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", names (l_gr_annotation_tr_bg[[i]][j]))
#     d_tr <- DataTrack(GR, name = id, background.title = cb_palette[i],
#                       type="heatmap", ylim = c(0, 0.5),
#                       gradient=c('white','blue'))#, fill=col_ctrl, background.title = col_ctrl) 
# 
#     list_all_bg <- append (list_all_bg, d_tr)
#   }
#   
# }

l_all_common_int <- list() 
common_intervals <- Reduce(subsetByOverlaps, c(unlist (l_granges_bg))) 

l_gr_annotation_tr_bg <- l_granges_bg

group_lab <- unlist(lapply (seq_along(l_granges_bg), function (i_group_exp) {
                                  rep (names(l_granges_bg[i_group_exp]), length(l_granges_bg[[i_group_exp]])) 
             }))

# for (i in 1:length(l_gr_annotation_tr_bg)){  
#   l_gr_common_int <- sapply (l_gr_annotation_tr_bg[[i]], function (l, common_GR=common_intervals) { 
#     mcol <- mcols(subsetByOverlaps (l, common_intervals)) 
#     return (mcol)
# #     return (data.frame(mcol))
#   })  
# #   l_all_common_int <- cbind(l_all_common_int, l_gr_common_int)  
#   l_all_common_int <- c(l_all_common_int, l_gr_common_int)  
# }

l_all_common_int <- sapply(unlist(l_gr_annotation_tr_bg), 
                           function (l, common_GR=common_intervals) { 
                            mcol <- mcols(subsetByOverlaps (l, common_intervals)) 
                            return (mcol)
#                             return (data.frame(mcol))
})

df <- as.data.frame (unlist(l_all_common_int))

# This was not working problably because number of rows was not correctly set
# df <- data.frame(matrix(unlist(l_all_common_int), nrow=length(common_intervals), byrow=T))
# names(df) <- paste ("id_", gsub(".+tr_(\\d+)(_.+$)", "\\1", names (unlist(l_gr_annotation_tr_bg))), sep="")
id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", names (unlist(l_gr_annotation_tr_bg)))
data_type <-  gsub(".+tr_(\\d+)_dt_(\\w+._.+$)", "\\2", names (unlist(l_gr_annotation_tr_bg)))
names(df) <- paste ("id", id, data_type, sep="_")

gr_common_intervals <- GRanges()
gr_common_intervals <- common_intervals
mcols(gr_common_intervals) <- df

# gr_common_intervals[, which(group_lab=="control")]


#####
## Problem with order of colors, the colors are not set by the provided order by the
## alphabetical order of the groups label, for instance if we have control and case
## the case color will be the first on the col assignment 
group_lab <- factor(group_lab, levels = unique(group_lab))
color_by_tr <- unlist(l_gr_color[unique(group_lab)])

# common_bedg_dt <- DataTrack(gr_common_intervals, name = "mean intake (mg)", type = "a",
#                                     showSampleNames = TRUE, #ylim = c(0, 0.5),                                     
#                                     groups = group_lab, col = color_by_tr,
#                                     background.title = col_back_title, size = tr_sum_size,
#                                     legend = TRUE)
# plotTracks(common_bedg_dt)#comment

# common_bedg_dt_boxPlot <- DataTrack(gr_common_intervals, name = "mean intake (mg)", type="a",
#                   showSampleNames = TRUE, #ylim = c(0, 0.5),                                     
#                   groups = group_lab, col=color_by_tr,
#                   legend=FALSE)

g_tr <- GenomeAxisTrack()

shinyServer(function(input, output) {
#   output$genomicPositionSelect <- renderUI({
    #     sliderInput( "tpos", "Time Point:", min = 10, max = g_max_end - 10, value = g_min_start + 10 )
#     sliderInput( "tpos", "Time Point:", min = 0, max = g_max_end - 10, value = 0 )
#   })
  
  #   pos <-  reactive({
  #     min( max( input$windowsize + 1, input$tpos ), max(g_max_end) - input$windowsize - 1 )    
  #   })
  
#   output$windowsize <- renderUI({                                                             
#     sliderInput("windowsize", "Window size:", min = min(g_min_start, 1000), max = max(g_max_end, 1000000), 
#                 value =min(g_max_end, 3000), step = min(g_max_end, 300))
#   })
  output$dataInterval <- renderUI({
    sliderInput("dataInterval", "Data interval:", 
                min = min(g_min_start, 1000), max = max(g_max_end, 1000000), 
                value = c(min(g_min_start, 1000), g_min_start+5000), step= 1000)
  }) 
  output$bedGraphRange <- renderUI({
    sliderInput("bedGraphRange", "Range bedgraph:", 
                min = min_v, max = max_v, value = c(0, 0.5), step= 0.1)
  })
  output$groups <- renderUI({
    checkboxGroupInput( "groups", "Groups to render", choices = unique(group_lab), selected=unique(group_lab))
  })
  output$groups_plot <- renderUI({                                                             
    checkboxInput("groups_plot", "Add group plot", FALSE)
  })
  output$boxplot <- renderUI({                                                             
    checkboxInput("boxplot", "Add boxplot:", FALSE)
  })  
  
  groups_dt <- reactive({
    if(!is.null(input$groups_plot) && input$groups_plot == TRUE) {
      
      gr_common_intervals_subset <- gr_common_intervals [ , which(group_lab==input$groups)]
      
      common_bedg_dt <- DataTrack(gr_common_intervals_subset, name = "mean intake (mg)", type = "a",
                showSampleNames = TRUE, #ylim = c(0, 0.5),                                     
                groups = group_lab[which(group_lab==input$groups)], col = color_by_tr,
                background.title = col_back_title, size = tr_sum_size,
                legend = TRUE)
      common_bedg_dt
    }
  })
  
  #  boxplot datatrack
  # it is overlap and then is very difficult to see anything
  boxplot_dt <- reactive({
    if(!is.null(input$boxplot) && input$boxplot == TRUE) {
      common_bedg_dt_boxplot <- groups_dt()
      displayPars(common_bedg_dt_boxplot) <- list(type=c("boxplot"))
      
      #       for (i in 1:length(list_gr)){
      #         displayPars(list_gr[[i]]) <- list(type=c("boxplot"), fill=cb_palette[i])
      # #           list(type=c("boxplot"), fill=cb_palette[i])
      #       }
      #       
      #       o_tr_boxplot <-OverlayTrack(list_gr)
      #       
      #       IdeogramTrack(genome=input$ucscgen, chromosome=input$chr,
      #                     showId=TRUE, showBandId=TRUE)
      #       o_tr_boxplot
      common_bedg_dt_boxplot
    }
  })
  
#   output$text1 <- renderText({ 
# #     paste (as.character (input$boxplot))
#     paste (as.character (input$groups)) 
#   })
  
  output$plotbed <- renderPlot({
#     if(length(input$windowsize)==0){
    if(length(input$bedGraphRange)==0){
      return(NULL)
    }
    else{
      if (input$boxplot==FALSE){
#         pt <- plotTracks(c(g_tr, list_all, list_all_bg, common_bedg_dt), 
#         pt <- plotTracks(c(g_tr, list_all, list_all_bg, groups_dt()), 
        pt <- plotTracks(c(g_tr, unlist(l_gr_annotation_tr_bed[input$groups]), 
                           unlist(list_all_bg[input$groups]), groups_dt()),
#         pt <- plotTracks(c(g_tr, unlist(l_gr_annotation_tr_bed[c("case", "control")]), list_all_bg, groups_dt()),
#         pt <- plotTracks(c(g_tr, list_all, o_tr),
#                          from=pos(), to=pos() + input$windowsize,
#                          from=input$tpos, to=input$tpos+ input$windowsize,
                         from=input$dataInterval[1], to=input$dataInterval[2], 
                         ylim=c(input$bedGraphRange[1], input$bedGraphRange[2]),                                                      
                         shape = "box", stacking = "dense")        
      }
      else {
        
#         pt <- plotTracks(c(g_tr, list_all, list_all_bg, common_bedg_dt, boxplot_dt()), 
#         pt <- plotTracks(c(g_tr, list_all, list_all_bg, groups_dt(), boxplot_dt()),
        pt <- plotTracks(c(g_tr, unlist(l_gr_annotation_tr_bed[input$groups]), unlist(list_all_bg[input$groups]), groups_dt(), boxplot_dt()),
#       pt <- plotTracks(c(g_tr, list_all, o_tr),
#                          from=pos(), to=pos() + input$windowsize,
#                          from=input$tpos, to=input$tpos+ input$windowsize,
                         from=input$dataInterval[1], to=input$dataInterval[2],
                         ylim=c(input$bedGraphRange[1], input$bedGraphRange[2]),
                         shape = "box", stacking = "dense")
      }
      
      return(pt)
    }
  })
})