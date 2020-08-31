params <- c(
  "Net turnover",
  "Extinction fraction",
  "Speciation",
  "Extinction",
  "Net diversification"
)

theme_h_scatter <- theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.background = element_rect(color="black", size = 1))

##### --- m_scatterplot_cp ui  ---------------------- #####

m_scatterplot_cp_ui <- function(id) {
  ns <- NS(id)
  tagList(
    checkboxGroupButtons(
      inputId = ns("m_scatter_cp"),
      choiceNames = "Scatterplot of model averaged diversification rates with two-dimensional colorplane",
      choiceValues = 1,
      status = "primary",
      selected = "Scatterplot of model averaged diversification rates with two-dimensional colorplane"
    ),
    
    conditionalPanel(
      condition = paste0("input['", ns("m_scatter_cp"), "'] == 1"),
      wellPanel(fluidRow(
        column(
          3,
          # parameter
          selectInput(
            inputId = ns("parameter"),
            label = "Parameter:",
            choices = params,
            selected = "Net turnover"
          ),
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
            label = "Secon character name:",
            placeholder = "p(x0)",
            value = "p(x0)"
          ),
          # plot_as_waiting_time
          checkboxInput(
            inputId = ns("plot_as_waiting_time"),
            label = "Plot as waiting time",
            value = FALSE
          ),
          actionButton(inputId = ns("plot"), label = "Plot"),
          actionButton(inputId = ns("code"), label = "Code")
        ),
        column(width = 4, 
               plotOutput(ns("plt")),
               tags$hr(),
               uiOutput(ns("plt_txt"), container = tags$code)
        )
      )
      )))
}

##### --- m_scatterplot_cp srv ---------------------- #####

m_scatterplot_cp_srv <-
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
      p <- m_scatterplot_cp(
        processed_recon = h_proc(),
        parameter = param(),
        focal_character = input$focal_char,
        focal_character_label = input$fc_label,
        second_character_label = input$sc_label,
        plot_as_waiting_time = input$plot_as_waiting_time,
        colors= c("#21908CFF", "#440154FF", "#FDE725FF")
      ) +
        theme_h_scatter
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
                         "<br/>m_scatterplot_cp(",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;processed_recon = m_proc,",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;focal_character = '", input$focal_char,"',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;focal_character_label = '", input$fc_label, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;second_character_label = '", input$sc_label, "',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;colors = c('#21908CFF', '#440154FF', '#FDE725FF'),", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;plot_as_waiting_time = ", input$plot_as_waiting_time,
                         "<br/>) + theme_classic()",
                         "<br/>",
                         "<br/># For more information see ?utilhisse::m_scatterplot_cp", sep="")
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