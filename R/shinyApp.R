library(shiny)
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Plot Generated Fuel Maps"),
  
  # Sidebar layout with input and output definitions ----
  
  fluidRow(
    column(2,
           numericInput("dimX","X dimention (m)",30,1,500,1)
    ),
    column(2, offset = 0,
           numericInput("dimY","Y dimention (m)",30,1,500,1)
    )
  ),
  
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "rho",
                  label = "heterogeneity:",
                  min = .1,
                  max = 10,
                  value = 5, step=.1),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "heterogeneity_scale",
                  label = "heterogeneity scale:",
                  min = .01,
                  max = 1,
                  value = .5, step=.01),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "density",
                  label = "density (fuel/m^2):",
                  min = .01,
                  max = 1,
                  value = .2, step=.01),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "radius",
                  label = "radius mean (m):",
                  min = .01,
                  max = 3,
                  value = 1, step=.1),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "sd_radius",
                  label = "radius stdev (m):",
                  min = 0,
                  max = 1,
                  value = .25, step=.05),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "height",
                  label = "height mean (m):",
                  min = .01,
                  max = 3,
                  value = 1, step=.1),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "sd_height",
                  label = "height stdev (m):",
                  min = 0,
                  max = 1,
                  value = .1, step=.05),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "reps",
                  label = "replicates:",
                  min = 1,
                  max = 25,
                  value = 9, step=1, round=T),
      actionButton("runcode", "Run"),
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "fuelsplot")
      
    )
  ),
  
  fluidRow(
    column(1,
           actionButton("loadchm", "Load CHM")
    ),
    column(2,
           textInput("chmfile",label=NULL,value = here::here("examples/CHM/LA_CHM.tif"), placeholder="path/to/canopy_height_model.tif")
    ),
    column(1,
           actionButton("rmchm", "Remove CHM")
    )
  ),
  fluidRow(
    column(2,
           numericInput("CHM_B","CHM scale [-1,1]",value=0,min=-5,max=5)
    )
  ),
  
  fluidRow(
    column(1,
           actionButton("expdatabutton", "Export Data")
    ),
    column(4, offset = 0,
           textInput("expdatafile", NULL, placeholder="path/to/datafile.csv")
    )
  ),
  
  fluidRow(
    column(1,
           actionButton("expplotbutton", "Export Plot")
    ),
    column(4, offset = 0,
           textInput("expplotfile", NULL, placeholder="path/to/plotfile.pdf")
    )
  )
)

server <- function(input, output, session) {
  
  chm <<- FALSE
  
  observeEvent(input$runcode, {
    if(chm){
      cat('Generating fuels with Canopy Height Model, B =',input$CHM_B,'.\n')
      fuel <<- fuelsgen::gen_fuels(dimX=isolate(input$dimX),dimY=isolate(input$dimY),
                                   density=isolate(input$density),heterogeneity=isolate(input$rho),
                                   radius=isolate(input$radius),sd_radius=isolate(input$sd_radius),
                                   height=isolate(input$height),sd_height=isolate(input$sd_height),
                                   reps = isolate(input$reps),
                                   sp.cov.locs = list(chm_locs),
                                   sp.cov.vals = list(chm_heights),
                                   sp.cov.scale = isolate(matrix(input$CHM_B)),
                                   heterogeneity.scale = isolate(input$heterogeneity_scale),verbose=T)
    } else{
      fuel <<- fuelsgen::gen_fuels(dimX=isolate(input$dimX),dimY=isolate(input$dimY),
                                   density=isolate(input$density),heterogeneity=isolate(input$rho),
                                   radius=isolate(input$radius),sd_radius=isolate(input$sd_radius),
                                   height=isolate(input$height),sd_height=isolate(input$sd_height),
                                   reps = isolate(input$reps),
                                   heterogeneity.scale = isolate(input$heterogeneity_scale),verbose=T)
    }
    
    output$fuelsplot <- renderPlot({
      plots <<- fuelsgen::plot_fuels(fuel)
      print(plots)
    }, height=1000,width=1000)
  })
  
  observeEvent(input$loadchm, {
    if(input$chmfile != ""){
      # load chm
      chm = raster::raster(input$chmfile)
      # overwrite user dimX, dimY when reading in CHM
      updateNumericInput(session,"dimX",value=chm@extent@xmax-chm@extent@xmin)
      updateNumericInput(session,"dimY",value=chm@extent@ymax-chm@extent@ymin)
      raster::extent(chm) = c(0,chm@extent@xmax-chm@extent@xmin,0,chm@extent@ymax-chm@extent@ymin)
      coords = raster::coordinates(chm)
      sp.cov.locs = list()
      sp.cov.locs$x = unique(coords[,1])
      sp.cov.locs$y = rev(unique(coords[,2])) # y's must be in increasing order for pracma::interp2
      canopy_heights = raster::values(chm,format='matrix')
      canopy_heights[canopy_heights<0] = 0
      canopy_heights = (canopy_heights - min(canopy_heights)) / (max(canopy_heights) - min(canopy_heights))
      canopy_heights = 2*canopy_heights - 1 # scale to [-1,1] so that low canopy reduces probability of shrub
      # flip y dim to match sp.cov.locs$y
      canopy_heights = canopy_heights[nrow(canopy_heights):1,]
      chm_heights <<- canopy_heights
      chm_locs <<- sp.cov.locs
      chm <<- TRUE
    }
  })
  observeEvent(input$rmchm, {
      chm <<- FALSE
  })
  observeEvent(input$expdatabutton, {
    if(input$expdatafile != "")
      # write fuels to csv
      if(input$reps>1){
        df = data.frame()
        for(i in 1:length(fuel$dat)){
          fuel$dat[[i]]$rep = i
          df = rbind(df,fuel$dat[[i]])
        }
        write.csv(df, file=input$expdatafile, row.names = F)
      } else{
        write.csv(fuel$dat, file=input$expdatafile, row.names = F)
      }
  })
  
  observeEvent(input$expplotbutton, {
    if(input$expplotfile != ""){
      pdf(input$expplotfile)
      print(plots)
      dev.off()
    }
  })
  
}

shinyApp(ui = ui, server = server)
