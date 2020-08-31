params <- c(
  "Net turnover",
  "Extinction fraction",
  "Speciation",
  "Extinction",
  "Net diversification"
)

theme_h_dot <- theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.background = element_rect(color="black", size = 1))

##### --- m_dotplot ui  ---------------------- #####

m_dotplot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    checkboxGroupButtons(
      inputId = ns("m_dot"),
      choiceNames = "Stacked dotplot of model averaged diversification rates",
      choiceValues = 1,
      status = "primary",
      selected = "Stacked dotplot of model averaged diversification rates"
    ),
    
    conditionalPanel(
      condition = paste0("input['", ns("m_dot"), "'] == 1"),
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
          textInput(
            inputId = ns("states_names1"),
            label = "State 00:",
            placeholder = "Character state 00",
            value = "00"
          ),
          textInput(
            inputId = ns("states_names2"),
            label = "State 01:",
            placeholder = "Character state 01",
            value = "01"
          ),
          textInput(
            inputId = ns("states_names3"),
            label = "State 10:",
            placeholder = "Character state 10",
            value = "10"
          ),
          textInput(
            inputId = ns("states_names4"),
            label = "State 11:",
            placeholder = "Character state 11",
            value = "11"
          ),
          numericInput(
            inputId = ns("bin_width"),
            label = "Bin width:",
            value = 0.05, 
            min=0,
            step = 0.01
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
        column(width = 9, 
               plotOutput(ns("plt")),
               tags$hr(),
               uiOutput(ns("plt_txt"), container = tags$code)
        )
      )
      )))
}

##### --- m_dotplot srv ---------------------- #####

m_dotplot_srv <-
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
      p <- m_dotplot(
        processed_recon = h_proc(),
        parameter = param(),
        states_names = c(input$states_names1, input$states_names2, input$states_names3, input$states_names4),
        plot_as_waiting_time = input$plot_as_waiting_time,
        bin_width = input$bin_width,
        colors= viridis(end = 0.6, n=4)
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
                         "<br/>h_proc <- h_process_recon(your_hisse_recon_object)",
                         "<br/>h_dotplot(",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;processed_recon = h_proc,",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;parameter = '", c(param()),"',",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;states_names = 
                         c('", input$states_names1, "','", input$states_names2, "','", input$states_names3, "','", input$states_names4, "'),",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;bin_width = ", input$bin_width, ",",
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;colors = c('#440154FF', '#414487FF', '#2A788EFF', '#22A884FF'),", 
                         "<br/>&nbsp;&nbsp;&nbsp;&nbsp;plot_as_waiting_time = ", input$plot_as_waiting_time,
                         "<br/>) + theme_classic()",
                         "<br/>",
                         "<br/># For more information see ?utilhisse::m_dotplot", sep="")
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