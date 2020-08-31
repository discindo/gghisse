##### --- h_trait_recon ui  ---------------------- #####

h_trait_recon_ui <- function(id) {
  ns <- NS(id)
  tagList(
    checkboxGroupButtons(
      inputId = ns("h_trait"),
      choiceNames = "Tree plot with ancestral reconstruction for the character states",
      choiceValues = 1,
      status = "primary",
      selected = "Tree plot with ancestral reconstruction for the character states"
    ),
    
    conditionalPanel(
      condition = paste0("input['", ns("h_trait"), "'] == 1"),
      wellPanel(fluidRow(
      column(
        width = 3,
        # x_label
        textInput(
          inputId = ns("x_label"),
          label = "Trait name:",
          placeholder = "The binary trait in your model"
        ),
        textInput(
          inputId = ns("states_names1"),
          label = "State 0:",
          placeholder = "Character state 0",
          value = 0
        ),
        textInput(
          inputId = ns("states_names2"),
          label = "State 1:",
          placeholder = "Character state 1",
          value = 1
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
          numericInput(
            inputId = ns("cutoff"),
            label = "Cut point in the range 0-1",
            min = 0,
            max = 1,
            value = .5,
            step = 0.1
          )
        ),
        actionButton(inputId = ns("plot"), label = "Plot"),
        actionButton(inputId = ns("code"), label = "Code")
      ),
      column(9,
             plotOutput(ns("plt"), height = 1000),
             tags$hr(),
             uiOutput(ns("plt_txt"), container = tags$code))
    )
  )))
}

##### --- h_trait_recon srv ---------------------- #####

h_trait_recon_srv <-
  function(input,
           output,
           session,
           h_obj) {
    h_proc <- reactive({
      x <- h_process_recon(h_obj ())
      return(x)
    })
    
    plt <- eventReactive(input$plot, {
      p <- h_trait_recon(
        processed_recon = h_proc(),
        trait_name = input$x_label,
        show_tip_labels = input$show_tip_labels,
        discrete = input$discrete,
        cutoff = input$cutoff,
        tree_layout = input$tree_layout,
        tree_direction = input$tree_direction,
        time_axis_ticks = input$time_axis_ticks,
        open_angle = input$open_angle,
        colors = c("firebrick", "steelblue")
      ) + theme(plot.background = element_rect(color="black", size = 1))
      return(p)
    })
    
    output$plt <- renderPlot({
      plt()
    })
    
    plt_txt <- eventReactive(input$code, {
      
      code_text <- paste("<b>Code to reproduce this figure in an R session: </b><br/>",
                         "<br/>",
                         "library(hisse)",
                         "<br/>library(utilhisse) # will load other dependencies",
                         "<br/>h_proc <- h_process_recon(your_hisse_recon_object)",
                         
                         "<br/>h_trait_recon(",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;processed_recon = h_proc,",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;trait_name = '", input$x_label, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;states_names = c(", input$states_names1, ",", input$states_names2, "),",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;show_tip_labels = ", input$show_tip_labels, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;discrete = ", input$discrete, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;cutoff = ", input$cutoff, ",", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;tree_layout = '", input$tree_layout, "',", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;tree_direction = '", input$tree_direction, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;time_axis_ticks = ", input$time_axis_ticks, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;open_angle = ", input$open_angle, ",", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;colors = c('firebrick', 'steelblue')",
                         "<br/>)",
                         "<br/># For more information see ?utilhisse::h_trait_recon", sep="")
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