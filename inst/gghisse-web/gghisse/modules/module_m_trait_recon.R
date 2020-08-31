# cols <- c("#00204DFF", "#1B3B6DFF", "#4E576CFF", "#727274FF", "#958F78FF", "#BCAF6FFF", "#E7D159FF")
cols <- viridis(n=7)
cols <- sample(cols, size = 7, replace = FALSE) 

##### --- m_trait_recon ui  ---------------------- #####

m_trait_recon_ui <- function(id) {
  ns <- NS(id)
  tagList(
    checkboxGroupButtons(
      inputId = ns("m_trait"),
      choiceNames = "Tree plot with ancestral reconstruction for the character states",
      choiceValues = 1,
      status = "primary",
      selected = "Tree plot with ancestral reconstruction for the character states"
    ),
    
    conditionalPanel(
      condition = paste0("input['", ns("m_trait"), "'] == 1"),
      wellPanel(fluidRow(
        column(
          width = 3,
          # first char states
          textInput(
            inputId = ns("fc_states1"),
            label = "States of the first character:",
            placeholder = "Name of state 0",
            value = "0"
          ),
          textInput(
            inputId = ns("fc_states2"),
            label = "",
            placeholder = "Name of state 1",
            value = "1"
          ),
          # second char states
          textInput(
            inputId = ns("sc_states1"),
            label = "States of the second character:",
            placeholder = "Name of state 0",
            value = "0"
          ),
          textInput(
            inputId = ns("sc_states2"),
            label = "",
            placeholder = "Name of state 1",
            value = "1"
          ),
          numericInput(
            inputId = ns("fc_cutoff"),
            label = "Cutoff to discretize probabilities for the first character:",
            min = 0,
            max = 1,
            value = 0.5,
            step = 0.1
          ),
          numericInput(
            inputId = ns("sc_cutoff"),
            label = "Cutoff to discretize probabilities for the second character:",
            min = 0,
            max = 1,
            value = 0.5,
            step = 0.1
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
          actionButton(inputId = ns("plot"), label = "Plot"),
          actionButton(inputId = ns("code"), label = "Code")
        ),
        column(width = 9, 
               plotOutput(ns("plt"), height = 1000),
               tags$hr(),
               uiOutput(ns("plt_txt"), container = tags$code)
        )
      )
      )))
}

##### --- n_trait_recon srv ---------------------- #####

m_trait_recon_srv <-
  function(input,
           output,
           session,
           h_obj) {
    h_proc <- reactive({
      x <- m_process_recon(h_obj ())
      return(x)
    })
    
    plt <- eventReactive(input$plot, {
      p <- m_trait_recon(
        processed_recon = h_proc(),
        states_of_first_character = c(input$fc_states1, input$fc_states2),
        states_of_second_character = c(input$sc_states1, input$sc_states2),
        show_tip_labels = input$show_tip_labels,
        cutoff = as.numeric(c(input$fc_cutoff, input$sc_cutoff)),
        tree_layout = input$tree_layout,
        tree_direction = input$tree_direction,
        time_axis_ticks = input$time_axis_ticks,
        open_angle = input$open_angle,
        colors= cols
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
                         "<br/>m_proc <- m_process_recon(your_hisse_recon_object)",
                         "<br/>m_trait_recon(",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;processed_recon = m_proc,",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;states_of_first_character = c('", input$fc_states1, "','", input$fc_states2, "'),",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;states_of_second_character = c('", input$sc_states1, "','", input$sc_states2, "'),",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;cutoff = as.numeric(c('", input$fc_cutoff, "','", input$sc_cutoff, "')),",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;colors = c('#440154FF', '#8FD744FF', '#21908CFF', '#443A83FF', '#FDE725FF', '#31688EFF', '#35B779FF'),",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;show_tip_labels = ", input$show_tip_labels, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;tree_layout = '", input$tree_layout, "',", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;tree_direction = '", input$tree_direction, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;time_axis_ticks = ", input$time_axis_ticks, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;open_angle = ", input$open_angle,
                         "<br/>)",
                         "<br/>",
                         "<br/># For more information see ?utilhisse::m_trait_recon", sep="")
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