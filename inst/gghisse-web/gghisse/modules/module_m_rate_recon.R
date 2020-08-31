params <- c(
  "Net turnover",
  "Extinction fraction",
  "Speciation",
  "Extinction",
  "Net diversification"
)

##### --- m_rate_recon ui  ---------------------- #####

m_rate_recon_ui <- function(id) {
  ns <- NS(id)
  tagList(
    checkboxGroupButtons(
      inputId = ns("m_rate"),
      choiceNames = "Tree plot with ancestral reconstruction for diversification rates",
      choiceValues = 1,
      status = "primary",
      selected = "Tree plot with ancestral reconstruction for diversification rates"
    ),
    
    conditionalPanel(
      condition = paste0("input['", ns("m_rate"), "'] == 1"),
      wellPanel(fluidRow(
        column(
          width = 3,
          # parameter
          selectInput(
            inputId = ns("parameter"),
            label = "Parameter:",
            choices = params,
            selected = "Net turnover"
          ),
          selectInput(
            inputId = ns("tree_layout"),
            label = "Tree layout:",
            choices = c('rectangular', 'slanted', 'circular', 'fan', 'radial'),
            selected = "slanted"
          ),
          selectInput(
            inputId = ns("tree_direction"),
            label = "Tree direction:",
            choices = c('right', 'left', 'up', 'down'),
            selected = "right"
          ),
          # show tip labels
          checkboxInput(
            inputId = ns("show_tip_labels"),
            label = "Show tip labels (illegible for large trees)",
            value = FALSE
          ),
          numericInput(
            inputId = ns("time_axis_ticks"),
            label = "Number of ticks for the time axis:",
            min = 0,
            max = 20,
            value = 10,
            step = 1
          ),
          numericInput(
            inputId = ns("open_angle"),
            label = "Space in degrees between the first and last tip (when tree layout is 'fan'):",
            min = 0,
            max = 360,
            value = 10,
            step = 5
          ),
          # discrete
          checkboxInput(
            inputId = ns("discrete"),
            label = "Discretize the probabilities",
            value = FALSE
          ),
          # cutoff
          conditionalPanel(
            condition = paste0("input['", ns("discrete"), "']"),
            helpText("Set a sequence of break points:"),
            numericInput(
              inputId = ns("begin"),
              label = "Min:",
              min = 0,
              value = 0,
              step = 0.1
            ),
            numericInput(
              inputId = ns("end"),
              label = "Max:",
              min = 0,
              value = 2,
              step = 0.1
            ),
            numericInput(
              inputId = ns("step"),
              label = "Step:",
              min = 0.01,
              value = 0.3,
              step = 0.1
            )
          ),
          actionButton(inputId = ns("plot"), label = "Plot"),
          actionButton(inputId = ns("code"), label = "Code")
        ),
        column(width = 9, 
               plotOutput(ns("plt"), height = 1000),
               tags$hr(),
               uiOutput(ns("plt_txt"), container = tags$code))
      )
      )))
}

##### --- m_rate_recon srv ---------------------- #####

m_rate_recon_srv <-
  function(input,
           output,
           session,
           h_obj) {
    h_proc <- reactive({
      x <- m_process_recon(h_obj ())
      return(x)
    })
    
    param <- reactive({
      x <- case_when(
        input$parameter == "Net turnover" ~ "turnover",
        input$parameter == "Extinction fraction" ~ "extinct.frac",
        input$parameter == "Speciation" ~ "speciation",
        input$parameter == "Extinction" ~ "extinction",
        input$parameter == "Net diversification" ~ "net.div"
      )
      return(x)
    })
    
    plt <- eventReactive(input$plot, {
      if (!input$discrete) {
        brk <- 1:2
      } else {
        brk <- seq(input$begin, input$end, input$step)
      }
      
      p <- m_rate_recon(
        processed_recon = h_proc(),
        parameter = param(),
        show_tip_labels = input$show_tip_labels,
        discrete = input$discrete,
        breaks =  seq(input$begin, input$end, input$step),
        tree_layout = input$tree_layout,
        tree_direction = input$tree_direction,
        time_axis_ticks = input$time_axis_ticks,
        open_angle = input$open_angle,
        colors = viridis(n=length(brk))
      ) + theme(plot.background = element_rect(color="black", size = 1))
      return(p)
    })
    
    output$plt <- renderPlot({
      plt()
    })
    
    
    plt_txt <- eventReactive(input$code, {
      if (!input$discrete) {
        brk <- 1:2
      } else {
        brk <- seq(input$begin, input$end, input$step)
      }
      code_text <- paste("<b>Code to reproduce this figure in an R session: </b><br/>",
                         "<br/>",
                         "library(hisse)",
                         "<br/>library(utilhisse) # will load other dependencies",
                         "<br/>m_proc <- m_process_recon(your_hisse_recon_object)",
                         
                         "<br/>m_rate_recon(",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;processed_recon = h_proc,",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;parameter = '", param(),"',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;show_tip_labels = ", input$show_tip_labels, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;discrete = ", input$discrete, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;breaks = seq(", input$begin, ",", input$end, ",", input$step, "),",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;plot_as_waiting_time = ", input$plot_as_waiting_time, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;tree_layout = '", input$tree_layout, "',", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;tree_direction = '", input$tree_direction, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;time_axis_ticks = ", input$time_axis_ticks, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;open_angle = ", input$open_angle, ",", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;colors = viridis(n = ", length(brk), ")",
                         "<br/>)",
                         "<br/># Note that `open_angle` and `breaks` are ignored unless `tree_layout = 'fan'` and `discrete=TRUE`, respectively",
                         "<br/># For more information see ?utilhisse::m_rate_recon", sep="")
      p <-
        wellPanel(
          class = "code_well",
          tags$style(".code_well {background-color: white ;}"),
          HTML(code_text)
        )
      return(p)
    })
    
    output$plt_txt <- renderUI({
      plt_txt()
    })
    
  }