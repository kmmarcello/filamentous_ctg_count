library(shiny)
library(shinydashboard)
library(shiny)
ui <- dashboardPage(
  dashboardHeader( title = "CTG count for Cyanobacteria Species"),
  dashboardSidebar(disable = T),
  dashboardBody(
      box(title = "Plot Options", width = 3,
          selectInput("x", "Select Gene", 
                      choices = c("gene_gyr_a", "gene_nus_a", "gene_inf_c","gene_inf_a", "gene_ots_a",
                                  "gene_dna_k", "gene_rec_a", "gene_dna_j","gene_ace_f", "gene_dea_d",
                                  "gene_inf_b", "gene_tig", "gene_rnr", "gene_dna_a", "gene_hup_b",
                                  "gene_rbf_a", "gene_yfl_a", "gene_pnp", "gene_csp", "gene_ace_e", "gene_des_a" ), 
                      selected= "gene_gyr_a")
    ),
    box(
      plotOutput("plot", width = "600px", height = "500px")
     )
  )
)

server <- function(input, output, session) {
  output$plot <- renderPlot({
    filament_cyanos %>% 
      filter(!is.na(tempurature_avg)) %>% 
      ggplot(aes_string(x = filament_cyanos$organism, y = input$x, fill = filament_cyanos$organism)) +
      geom_col()+
      theme(legend.position = "right",
            text = element_text(size=6),
            axis.text.x = element_text(angle = 60, hjust=1))+
      labs(title = "CTG count for Cyanobacteria Species",x = "Cyanobacteria", y = "Gene Count")
    
  })
  session$onSessionEnded(stopApp)
}
  


shinyApp(ui, server)







