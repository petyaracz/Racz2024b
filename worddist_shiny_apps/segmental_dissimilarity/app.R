# -- setup -- #

library(shiny)
library(dplyr)
library(readr)
library(knitr)
library(stringr)
library(glue)

# -- read -- #

# setwd('~/Github/packages/worddist/worddist_shiny_apps/segmental_dissimilarity/')

d_fm = read_tsv('files/app_hungarian.tsv')
d_nc = read_tsv('files/app_hungarian_nc.tsv')
d_dist = read_tsv('files/app_hungarian_nc_dist_w.tsv')

# -- fun -- #

# squeezeTable = function(dat){
#   dat |>
#     knitr::kable('simple') |>
#     capture.output() |>
#     paste(collapse = "\n")
# }

getDistance = function(s1, s2, d_dist){
  distance = d_dist |> 
    filter(segment1 == s1) |> 
    select(all_of(s2)) |> 
    pull() |> 
    round(2)
  
  glue('The distance is {distance} between /{s1}/ and /{s2}/.')
}


#   
#   non_shared_nc = all_nc |> 
#     anti_join(shared_nc, by = join_by(feature_bundle, segments))
#  

# -- ui -- #

# Define UI where you input two segments and see the output matrix and the calculated dissimilarity
ui <- fluidPage(

    # Application title
    titlePanel("Segmental dissimilarity calculator"),
    # Additional text
    p("- Based on Frisch 1997, Albright and Hayes 2003, Frisch, Pierrehumbert, and Broe 2004."),
    p("- Segmental dissimilarity is calculated as 1 - (the number of overlapping natural classes / (the number of overlapping and non-overlapping natural classes))."),
    p("- Natural classes are specified using minimal feature descriptions as per Albright 2003."),
    p("- The underspecified feature matrix is based on Siptár and Törkenczy 2000 with changes by Tóth and Rácz (in prep)."),
    p("Use the drop-down menus to pick the two segments which you want to calculate distance for. The right-hand panel is populated once you press the button."),
    p("Note: /c/ is IPA /ts/. /ç/ is IPA /c/. /ʝ/ is IPA /ɟ/. Hungarian /v/ is an underspecified obstruent."),

    # Sidebar with two rolldown menus for the two segments
    sidebarLayout(
        sidebarPanel(
            # First segment
            selectInput("segment1", "Segment 1:", choices = d_fm$segment),
            # Second segment
            selectInput("segment2", "Segment 2:", choices = d_fm$segment),
            actionButton("go", "Calculate Distance")
        ),
        
    
        # Main panel showing a table called overlaps
        mainPanel(
            h2("Segmental Dissimilarity Calculator"),
            h3("Dissimilarity"),
            textOutput("distance"),
            h4("Shared and non-shared classes"),
            p("Shared natural classes:"),
            tableOutput("nc_shared"),
            p("Non-shared natural classes:"),
            tableOutput("nc_not_shared"),
            h4("Natural classes for each"),
            textOutput("nc1_title"),
            tableOutput("nc1"),
            textOutput("nc2_title"),
            tableOutput("nc2"),
            h3("Additional Information"),
            # smaller text saying "Feature matrix"
            h4("Total feature matrix"),
            tableOutput("d_fm"),
            h4("Total distance matrix"),
            tableOutput("d_dist")
        )
    )
)

# -- server -- #

# Define server logic
server <- function(input, output) {

  observeEvent(input$go, {
    s1 = input$segment1
    s2 = input$segment2
    output$d_fm <- renderTable({ d_fm })
    output$d_dist <- renderTable({ d_dist })
    output$distance <- renderText({ getDistance(s1, s2, d_dist) })
    output$nc1_title <- renderText({ glue('Natural classes for /{s1}/:') })
    output$nc1 <- renderTable({filter(d_nc, str_detect(segments, s1))})
    output$nc2_title <- renderText({ glue('Natural classes for /{s2}/:') })
    output$nc2 <- renderTable({filter(d_nc, str_detect(segments, s2))})
    output$nc_shared <- renderTable({filter(d_nc, str_detect(segments, s1) & str_detect(segments, s2))})
    output$nc_not_shared <- renderTable({filter(d_nc, xor(str_detect(segments, s1), str_detect(segments, s2)))})
                              }
    )
  }

# Run the application 
shinyApp(ui = ui, server = server)
