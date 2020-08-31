library(shiny)
library(shinythemes)
library(shinyWidgets)
library(hisse)
library(gghisse)
library(viridis)
library(colorplaner)


refs_text <- tags$div(
  tags$h3("About"),
  tags$p(
    "This web application is envisioned as a companion to the R package `hisse` for Hidden Markov Model variants of State Speciation and Extinction models. The goal is to make it easier to quickly visualize HiSSE results and produce nice looking plots that with some customization can be used for reports or presentations. These functions do not replace hisse's plotting functions, rather provide alternative ways of looking at the data. All of the plotting functions used here are available from the package ",
    tags$a(href = "https://github.com/teofiln/utilhisse", 'utilhisse'),
    " and can be used indepently of this app."
  ),
  tags$p(
    "The workflow starts with fitting models in R using the functions `hisse` or `MuHiSSE` (GeoHiSSE models are not yet supported here) followed by reconstructing marginal ancestral states for traits and rates using the HiSSE functions `MarginRecon` and `MarginReconMuHiSSE`. The ancestral states objects (of class `hisse.states` or `muhisse.states`) can then be saved to an external file with R's `save` function. These saved reconstruction objects are then loaded into this app to visualize the results."
  ),
  tags$p(
    "This app was written by ",
    tags$a(href = "https://teofil.discindo.org", "Teofil Nakov"),
    ". The R Shiny code for the app is available on ",
    tags$a(href = "https://github.com/teofiln/hisse-web", "GitHub"),
    "."
  ),
  
  tags$h3("References"),
  tags$h4("Papers:"),
  tags$p(
    "Beaulieu, J. M., & O’Meara, B. C. (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic biology, 65(4), 583-601.",
    tags$a(href = "https://academic.oup.com/sysbio/article/65/4/583/1753616", "link")
  ),
  tags$p(
    "Caetano, D. S., O'Meara, B. C., & Beaulieu, J. M. (2018). Hidden state models improve state‐dependent diversification approaches, including biogeographical models. Evolution, 72(11), 2308-2324.",
    tags$a(href = "https://onlinelibrary.wiley.com/doi/abs/10.1111/evo.13602", "link")
  ),
  tags$p(
    "Nakov, T., Beaulieu, J. M., & Alverson, A. J. (2019). Diatoms diversify and turn over faster in freshwater than marine environments. bioRxiv, 406165.",
    tags$a(href = "https://doi.org/10.1101/406165", "link")
  ),
  tags$h4("R packages:"),
  tags$p(
    "Beaulieu, J., O'Meara, B., & Caetano, D. (2019). Package ‘hisse’",
    tags$a(href = "https://cran.r-project.org/web/packages/hisse/index.html", 'CRAN')
  ),
  tags$p(
    "Nakov, T. (2019). Package ’utilhisse’",
    tags$a(href = "https://github.com/teofiln/utilhisse", 'Guthub')
  )
)

# modules
mods <-
  list.files(path = "modules/",
             pattern = "*.R",
             full.names = TRUE)
lapply(mods, source)

ui <-
  shiny::navbarPage(
    theme = shinytheme("flatly"),
    title = "Hidden State Speciation and Extinction (HiSSE)",
    fluid = TRUE,
    inverse = TRUE,
    
    # navbar for type of model
    shiny::tabPanel(
      "BINARY",
      # file load
      shiny::column(3,
                    shiny::wellPanel(
                      shiny::fileInput(
                 width = "100%",
                 'h_recon_input',
                 'Choose a HiSSE marginal ancestral reconstruction file:',
                 accept = c('.Rsave', "RSave", '.Rdata', 'RData')
               ),
               shiny::helpText(
                 'This should be an object output from `hisse::MarginRecon` or a list of such objects where each element contains the AIC score of the model. Make sure the object was saved to a file with the extension "Rsave" or "Rdata".'
               ),
               shiny::checkboxInput("h_demo", label = "Use demo file", value = FALSE),
               shiny::conditionalPanel(condition = "input.h_demo", 
                                       shiny::HTML("This is an ancestral state reconstruction object for a CID4 model for plankton-benthos trait in diatoms. For more information see our paper: <a href=https://doi.org/10.1101/406165>'Biorxiv'</a>."))
             )),
      shiny::column(
        9,
        h_scatterplot_ui(id = "1"),
        h_dotplot_ui(id = "2"),
        h_ridgelines_ui(id = "3"),
        h_trait_recon_ui(id = "4"),
        h_rate_recon_ui(id = "5")
      )
      
    ),
    shiny::tabPanel(
      "MULTISTATE",
      # file load
      shiny::column(3,
                    shiny::wellPanel(
                      shiny::fileInput(
                 width = "100%",
                 'm_recon_input',
                 'Choose a MuHiSSE marginal ancestral reconstruction file:',
                 accept = c('.Rsave', "RSave", '.Rdata', 'RData')
               ),
               shiny::helpText(
                 'This should be an object output from `hisse::MarginReconMuHiSSE` or a list of such objects where each element contains the AIC score of the model. Make sure the object was saved to a file with the extension "Rsave" or "Rdata".'
               ),
               shiny::checkboxInput("m_demo", label = "Use demo file", value = FALSE),
               shiny::conditionalPanel(condition = "input.m_demo", 
                                HTML("This is an ancestral state reconstruction object for a MuHiSSE model for marine-freshwater + plankton-benthos interaction in diatoms. For more information see our paper: <a href=https://doi.org/10.1101/406165>'Biorxiv'</a>."))
    )),
    shiny::column(
        9,
        m_scatterplot_ui(id = "11"),
        m_dotplot_ui(id = "12"),
        m_ridgelines_ui(id = "13"),
        m_scatterplot_cp_ui(id = "14"),
        m_trait_recon_ui(id = "15"),
        m_trait_recon_cp_ui(id = "16"),
        m_rate_recon_ui(id = "17")
      )
    ),
    # tabPanel("GEOGRAPHIC")
    shiny::tabPanel("ABOUT",
                    shiny::column(width = 6,
                    wellPanel(refs_text)))
  )

server <- function(input, output) {
  ##### ---- server logic for binary ----------------------- #####
  h_recon_load <- shiny::reactive({
    if (input$h_demo) {
      data("diatoms")
      H <- diatoms$cid4_recon
      return(H)
    } else {
      shiny::validate(
        shiny::need(
          input$h_recon_input != "" ,
          "Please select a HiSSE ancestral reconstruction file"
        )
      )
      in_file <- input$h_recon_input
      
      H <- get(load(in_file$datapath))
      shiny::validate(
        shiny::need(
          class(H) == "hisse.states" || class(H[[1]]) == "hisse.states",
          "Looks like this is not a MuHiSSE ancestral reconstruction file (makes sure `class(obj)` or if list, `class(obj[[1]])`, returns 'hisse.states')"
        )
      )
      return(H)
    }
  })
  
  shiny::callModule(module = h_scatterplot_srv,
             id = "1",
             h_obj = h_recon_load)
  
  shiny::callModule(module = h_dotplot_srv,
             id = "2",
             h_obj = h_recon_load)
  
  shiny::callModule(module = h_ridgelines_srv,
             id = "3",
             h_obj = h_recon_load)
  
  shiny::callModule(module = h_trait_recon_srv,
             id = "4",
             h_obj = h_recon_load)
  
  shiny::callModule(module = h_rate_recon_srv,
             id = "5",
             h_obj = h_recon_load)
  
  
  ##### ---- server logic for multistate ----------------------- #####
  m_recon_load <- shiny::reactive({
    if (input$m_demo) {
      data("diatoms")
      H <- diatoms$muhisse_recon
      return(H)
    } else {
      shiny::validate(
        shiny::need(
          input$m_recon_input != "" ,
          "Please select a MuHiSSE ancestral reconstruction file"
        )
      )
      in_file <- input$m_recon_input
      H <- get(load(in_file$datapath))
      shiny::validate(
        shiny::need(
          class(H) == "muhisse.states" || class(H[[1]]) == "muhisse.states",
          "Looks like this is not a MuHiSSE ancestral reconstruction file (make sure `class(obj)` or if list, `class(obj[[1]])`, returns 'muhisse.states')"
        )
      )
      return(H)
    }
  })
  
  shiny::callModule(module = m_scatterplot_srv,
             id = "11",
             h_obj = m_recon_load)
  
  shiny::callModule(module = m_dotplot_srv,
             id = "12",
             h_obj = m_recon_load)
  
  shiny::callModule(module = m_ridgelines_srv,
             id = "13",
             h_obj = m_recon_load)
  
  shiny::callModule(module = m_scatterplot_cp_srv,
             id = "14",
             h_obj = m_recon_load)
  
  shiny::callModule(module = m_trait_recon_srv,
             id = "15",
             h_obj = m_recon_load)
  
  shiny::callModule(module = m_trait_recon_cp_srv,
             id = "16",
             h_obj = m_recon_load)

  shiny::callModule(module = m_rate_recon_srv,
             id = "17",
             h_obj = m_recon_load)
}

# Run the application
shiny::shinyApp(ui = ui, server = server)
