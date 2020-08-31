##### --- m_trait_recon_cp ui  ---------------------- #####

m_trait_recon_cp_ui <- function(id) {
  ns <- NS(id)
  tagList(
    checkboxGroupButtons(
      inputId = ns("m_trait_cp"),
      choiceNames = "Tree plot with ancestral reconstruction for the character states with two-dimensional colorplane",
      choiceValues = 1,
      status = "primary",
      selected = "Tree plot with ancestral reconstruction for the character states with two-dimensional colorplane"
    ),
    
    conditionalPanel(
      condition = paste0("input['", ns("m_trait_cp"), "'] == 1"),
      wellPanel(fluidRow(
        column(
          width = 3,
          # first char states
          selectInput(
            inputId = ns("focal_char"),
            label = "Focal character:",
            choices = c("prob_0x", "prob_x0"),
            selected = "prob_0x"
          ),
          textInput(
            inputId = ns("fc_label"),
            label = "Focal character name:",
            placeholder = "p(0x)",
            value = "p(0x)"
          ),
          textInput(
            inputId = ns("sc_label"),
            label = "Second character name:",
            placeholder = "p(x0)",
            value = "p(x0)"
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
          radioGroupButtons(
            inputId = ns("code"),
            label = "Code",
            choiceNames = "Code",
            choiceValues = 1
          ))
        
        ),
        column(width = 9, 
               plotOutput(ns("plt"), height = 1000),
               tags$hr(),
               uiOutput(ns("plt_txt"), container = tags$code)
        )
      )))
}

##### --- n_trait_recon_cp srv ---------------------- #####

m_trait_recon_cp_srv <-
  function(input,
           output,
           session,
           h_obj) {
    h_proc <- reactive({
      x <- m_process_recon(h_obj ())
      return(x)
    })
    
    plt <- eventReactive(input$plot, {
      p <- m_trait_recon_cp(
        processed_recon = h_proc(),
        focal_character = input$focal_char, 
        focal_character_label = input$fc_label,
        second_character_label = input$sc_label,
        show_tip_labels = input$show_tip_labels,
        tree_layout = input$tree_layout,
        tree_direction = input$tree_direction,
        time_axis_ticks = input$time_axis_ticks,
        open_angle = input$open_angle,
        colors= c("#21908CFF", "#440154FF", "#FDE725FF")
      ) + theme(plot.background = element_rect(color="black", size = 1))
      return(p)
    })
    
    output$plt <- renderPlot({
      plt()
    })
    
    
    plt_txt <- observeEvent(input$code == "Yes", {
      
      code_text <- paste("<b>Code to reproduce this figure in an R session: </b><br/>",
                         "<br/>",
                         "library(hisse)",
                          "<br/>library(utilhisse) # will load other dependencies",
                         "<br/>m_proc <- m_process_recon(your_hisse_recon_object)",
                         "<br/>m_trait_recon_cp(",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;processed_recon = m_proc,",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;focal_character = '", input$focal_char,"',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;focal_character_label = '", input$fc_label, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;second_character_label = '", input$sc_label, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;colors = c('#21908CFF', '#440154FF', '#FDE725FF'),",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;show_tip_labels = ", input$show_tip_labels, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;tree_layout = '", input$tree_layout, "',", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;tree_direction = '", input$tree_direction, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;time_axis_ticks = ", input$time_axis_ticks, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;open_angle = ", input$open_angle,
                         "<br/>)",
                         "<br/>",
                         "<br/># For more information see ?utilhisse::m_trait_recon_cp", sep="")
    
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